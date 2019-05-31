/*! Definition of a simple tuple producer.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
#include "h-tautau/Analysis/include/SimpleTuple.h"
#include "h-tautau/Instruments/include/SimpleTupleProducerConfigReader.h"
#include "h-tautau/McCorrections/include/TauIdWeight.h"

struct Arguments {
    REQ_ARG(std::string, cfg);
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
    OPT_ARG(unsigned, n_threads, 1);
};

namespace simple_tuple {

class SimpleTupleProducer {
public:
    SimpleTupleProducer(const Arguments& _args) :
        args(_args)
    {
        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());

        analysis::ConfigReader configReader;
        ProducerSetupCollection setups;
        ProducerSetupEntryReader setupReader(setups);
        configReader.AddEntryReader("SETUP", setupReader, false);

        SampleDescriptorEntryReader sampleReader(samples);
        configReader.AddEntryReader("SAMPLE", sampleReader, true);
        configReader.ReadConfig(args.cfg());

        if(setups.size() != 1)
            throw analysis::exception("To many setups in the config.");
        setup = setups.begin()->second;
    }

    void Run()
    {
        static constexpr float default_value = DefaultFillValue<float>();
        static const std::set<analysis::Channel> default_channels = {
            analysis::Channel::ETau, analysis::Channel::MuTau, analysis::Channel::TauTau
        };

        static const std::map<analysis::Channel, std::pair<double, double>> pt_cuts = {
            { analysis::Channel::ETau, { 26, 30 } },
            { analysis::Channel::MuTau, { 20, 30 } },
            { analysis::Channel::TauTau, { 50, 40 } },
        };


        // expected QCD/total_bkg in the signal region
        // eTau: 3.3 / 46.5 = 0.071
        // muTau: 5.0 / 77.7 = 0.064
        // tauTau: 35.5 / 156.2 = 0.228
        static const std::map<analysis::Channel, double> os_ss_sf = {
            { analysis::Channel::ETau, 1.10 }, //0.87 },
            { analysis::Channel::MuTau, 0.99 }, //0.7 },
            { analysis::Channel::TauTau, 1.20 }, //1.22 },
        };

        // static const std::map<analysis::Channel, double> os_ss_sf_in = {
        //     { analysis::Channel::ETau, 0.84 }, //0.87 },
        //     { analysis::Channel::MuTau, 0.70 }, //0.7 },
        //     { analysis::Channel::TauTau, 0.90 }, //1.22 },
        // };
        //
        // static const std::map<analysis::Channel, double> os_ss_sf_out = {
        //     { analysis::Channel::ETau, 1.18 }, //0.87 },
        //     { analysis::Channel::MuTau, 0.99 }, //0.7 },
        //     { analysis::Channel::TauTau, 1.27 }, //1.22 },
        // };


        static const double tau_id_sf = 0.87;
        analysis::mc_corrections::TauIdWeight tau_trig_weight_medium(
            "h-tautau/McCorrections/data/Tau/fitresults_tt_moriond2017.json", analysis::DiscriminatorWP::Medium);
        analysis::mc_corrections::TauIdWeight tau_trig_weight_tight(
            "h-tautau/McCorrections/data/Tau/fitresults_tt_moriond2017.json", analysis::DiscriminatorWP::Tight);

        auto output_file = root_ext::CreateRootFile(args.output());
        EventTuple tuple(output_file.get(), false);

        for(const auto& item : samples.get_ordered_by_insertion()) {
            const SampleDescriptor& sample = *item.second;

            std::cout << "Processing " << sample.name << "..." << std::endl;
            auto input_file = root_ext::OpenRootFile(args.input() + "/" + sample.file_path);
            auto summary_tuple = ntuple::CreateSummaryTuple("summary", input_file.get(), true,
                                                            ntuple::TreeState::Skimmed);
            const auto summary = ntuple::MergeSummaryTuple(*summary_tuple);

            const analysis::SummaryInfo summary_info(summary);

            std::set<analysis::Channel> channels = sample.channels;
            if(channels.empty())
                channels = default_channels;
            for(auto channel : channels) {
                auto input_tuple = ntuple::CreateEventTuple(ToString(channel), input_file.get(), true,
                                                            ntuple::TreeState::Skimmed);
                const auto& channel_triggers = setup.trigger_paths.at(channel);
                const auto& channel_pt_cuts = pt_cuts.at(channel);
                const double channel_os_ss_sf = os_ss_sf.at(channel);

                for(const auto& event : *input_tuple) {
                    if(!(event.p4_1.pt() > channel_pt_cuts.first && event.p4_2.pt() > channel_pt_cuts.second)) continue;

                    analysis::TriggerResults triggerResults;
                    triggerResults.SetAcceptBits(event.trigger_accepts);
                    triggerResults.SetMatchBits(event.trigger_matches);
                    triggerResults.SetDescriptors(summary_info.GetTriggerDescriptors(channel));
                    if(!triggerResults.AnyAcceptAndMatch(channel_triggers)) continue;

                    const bool os = event.q_1 * event.q_2 == -1;
                    if(!os && sample.sample_type == SampleType::H_tautau) continue;

                    double weight = 1;
                    if(sample.sample_type != SampleType::Data) {
                        const analysis::GenMatch gen_match_1 = static_cast<analysis::GenMatch>(event.gen_match_1);
                        const analysis::GenMatch gen_match_2 = static_cast<analysis::GenMatch>(event.gen_match_2);

                        weight = event.weight_total * sample.cross_section * setup.int_lumi
                               / summary.totalShapeWeight;

                        weight *= tau_trig_weight_tight.Get(event) / tau_trig_weight_medium.Get(event);
                        if(gen_match_1 == analysis::GenMatch::Tau)
                            weight *= tau_id_sf;
                        if(gen_match_2 == analysis::GenMatch::Tau)
                            weight *= tau_id_sf;
                    }
                    const double weight_sign = !os && sample.sample_type != SampleType::Data ? -1 : 1;
                    const double weight_sf = os ? 1 : channel_os_ss_sf;
                    weight = weight * weight_sign * weight_sf;
                    const SampleType sample_type = os ? sample.sample_type : SampleType::QCD;

                    tuple().event_id = GenerateEventIdentifier(event);
                    tuple().sample_type = static_cast<int>(sample_type);
                    tuple().higgs_prod_type = static_cast<int>(sample.higgs_prod);
                    tuple().higgs_gen_mass = sample.higgs_gen_mass;
                    tuple().channel_id = static_cast<int>(channel);
                    tuple().weight = static_cast<float>(weight);
                    tuple().npv = event.npv;
                    tuple().px_1 = event.p4_1.px();
                    tuple().py_1 = event.p4_1.py();
                    tuple().pz_1 = event.p4_1.pz();
                    tuple().E_1 = event.p4_1.E();
                    tuple().q_1 = event.q_1;
                    tuple().dxy_1 = event.dxy_1;
                    tuple().dz_1 = event.dz_1;
                    tuple().decay_mode_1 = event.decayMode_1;
                    tuple().gen_match_1 = event.gen_match_1;
                    tuple().gen_px_1 = HasGenMatch(event.gen_match_1) ? event.gen_p4_1.px() : default_value;
                    tuple().gen_py_1 = HasGenMatch(event.gen_match_1) ? event.gen_p4_1.py() : default_value;
                    tuple().gen_pz_1 = HasGenMatch(event.gen_match_1) ? event.gen_p4_1.pz() : default_value;
                    tuple().gen_E_1 = HasGenMatch(event.gen_match_1) ? event.gen_p4_1.E() : default_value;
                    tuple().px_2 = event.p4_2.px();
                    tuple().py_2 = event.p4_2.py();
                    tuple().pz_2 = event.p4_2.pz();
                    tuple().E_2 = event.p4_2.E();
                    tuple().q_2 = event.q_2;
                    tuple().dxy_2 = event.dxy_2;
                    tuple().dz_2 = event.dz_2;
                    tuple().decay_mode_2 = event.decayMode_2;
                    tuple().gen_match_2 = event.gen_match_2;
                    tuple().gen_px_2 = HasGenMatch(event.gen_match_2) ? event.gen_p4_2.px() : default_value;
                    tuple().gen_py_2 = HasGenMatch(event.gen_match_2) ? event.gen_p4_2.py() : default_value;
                    tuple().gen_pz_2 = HasGenMatch(event.gen_match_2) ? event.gen_p4_2.pz() : default_value;
                    tuple().gen_E_2 = HasGenMatch(event.gen_match_2) ? event.gen_p4_2.E() : default_value;
                    tuple().met_px = event.pfMET_p4.px();
                    tuple().met_py = event.pfMET_p4.py();
                    tuple().fitted_H_mass = event.SVfit_p4.mass();

                    auto jet_indices = OrderJetsByPt(event);
                    const size_t n_jets = event.jets_p4.size();
                    tuple().n_jets = n_jets;
                    tuple().leading_jet_px = n_jets > 0 ? event.jets_p4.at(jet_indices.at(0)).px() : default_value;
                    tuple().leading_jet_py = n_jets > 0 ? event.jets_p4.at(jet_indices.at(0)).py() : default_value;
                    tuple().leading_jet_pz = n_jets > 0 ? event.jets_p4.at(jet_indices.at(0)).pz() : default_value;
                    tuple().leading_jet_E = n_jets > 0 ? event.jets_p4.at(jet_indices.at(0)).E() : default_value;
                    tuple().leading_jet_btag = n_jets > 0 ? event.jets_csv.at(jet_indices.at(0)) : default_value;

                    tuple().subleading_jet_px = n_jets > 1 ? event.jets_p4.at(jet_indices.at(1)).px() : default_value;
                    tuple().subleading_jet_py = n_jets > 1 ? event.jets_p4.at(jet_indices.at(1)).py() : default_value;
                    tuple().subleading_jet_pz = n_jets > 1 ? event.jets_p4.at(jet_indices.at(1)).pz() : default_value;
                    tuple().subleading_jet_E = n_jets > 1 ? event.jets_p4.at(jet_indices.at(1)).E() : default_value;
                    tuple().subleading_jet_btag = n_jets > 1 ? event.jets_csv.at(jet_indices.at(1)) : default_value;

                    tuple.Fill();
                }
            }
        }

        tuple.Write();
    }

private:
    static ULong64_t GenerateEventIdentifier(const ntuple::Event& event)
    {
        size_t seed = std::hash<ULong64_t>{}(event.run);
        seed ^= std::hash<ULong64_t>{}(event.lumi) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= std::hash<ULong64_t>{}(event.evt) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
    }

    static bool HasGenMatch(Int_t gen_match)
    {
        return static_cast<analysis::GenMatch>(gen_match) != analysis::GenMatch::NoMatch;
    }

    static std::vector<size_t> OrderJetsByPt(const ntuple::Event& event)
    {
        std::vector<size_t> indices(event.jets_p4.size());
        std::iota(indices.begin(), indices.end(), 0);

        auto compare = [&](size_t a, size_t b) {
            return event.jets_p4.at(a).Pt() > event.jets_p4.at(b).Pt();
        };

        std::sort(indices.begin(), indices.end(), compare);
        return indices;
    }

private:
    Arguments args;
    ProducerSetup setup;
    SampleDescriptorCollection samples;
};

} // namespace simple_tuple

PROGRAM_MAIN(simple_tuple::SimpleTupleProducer, Arguments)
