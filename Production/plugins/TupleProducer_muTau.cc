/*! Implementation of an event tuple producer for the mu-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_muTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_muTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::MuTau;

    SelectionResults selection;
    cut(primaryVertex.isNonnull(), "vertex");

    if(applyTriggerMatch) cut(triggerTools.HaveTriggerFired(hltPaths), "trigger");

    // Di-Lepton Veto
    const auto z_muons = CollectZmuons();
    const auto z_muons_candidates = FindCompatibleObjects(z_muons, z_muons, ZmumuVeto::deltaR, "Z_mu_mu", 0);
    selection.Zveto = z_muons_candidates.size();

    // Signal-like leptons selection
    const auto selectedMuons = CollectSignalMuons();
    cut(selectedMuons.size(), "muons");
    const auto selectedTaus = CollectSignalTaus();
    cut(selectedTaus.size(), "taus");

    const double DeltaR_betweenSignalObjects = productionMode == ProductionMode::hh
            ? cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects
            : cuts::H_tautau_2016::DeltaR_betweenSignalObjects;
    auto higgses = FindCompatibleObjects(selectedMuons, selectedTaus, DeltaR_betweenSignalObjects, "H_mu_tau");
    cut(higgses.size(), "mu_tau_pair");

    std::vector<HiggsCandidate> selected_higgses = higgses;
    if(applyTriggerMatch) {
        selected_higgses = triggerTools.ApplyTriggerMatch(higgses, hltPaths, false);
        cut(selected_higgses.size(), "triggerMatch");
    }

    std::sort(selected_higgses.begin(), selected_higgses.end(), &HiggsComparitor<HiggsCandidate>);
    selection.SetHiggsCandidate(selected_higgses.front());

    //Third-Lepton Veto
    const auto electronVetoCollection = CollectVetoElectrons();
    const auto muonVetoCollection = CollectVetoMuons(&selection.higgs->GetFirstDaughter());
    selection.electronVeto = electronVetoCollection.size();
    selection.muonVeto = muonVetoCollection.size();

    ApplyBaseSelection(selection, selection.higgs->GetDaughterMomentums());
    if(runSVfit)
        selection.svfitResult = svfitProducer.Fit(*selection.higgs, *met);
    FillEventTuple(selection);
}

std::vector<BaseTupleProducer::MuonCandidate> TupleProducer_muTau::CollectSignalMuons()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_muTau::SelectSignalMuon, this, _1, _2);
    return CollectObjects("SignalMuons", base_selector, muons);
}

std::vector<BaseTupleProducer::TauCandidate> TupleProducer_muTau::CollectSignalTaus()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_muTau::SelectSignalTau, this, _1, _2);
    return CollectObjects("SignalTaus", base_selector, taus);
}

void TupleProducer_muTau::SelectSignalMuon(const MuonCandidate& muon, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::MuTau::muonID;

    cut(true, "gt0_mu_cand");
    const LorentzVector& p4 = muon.GetMomentum();
    double pt_cut = pt;
    if(productionMode == ProductionMode::hh) pt_cut = cuts::hh_bbtautau_2016::MuTau::muonID::pt;
    else if (productionMode == ProductionMode::h_tt_mssm) pt_cut = cuts::H_tautau_2016_mssm::MuTau::muonID::pt;
    cut(p4.pt() > pt_cut, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double muon_dxy = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muon_dxy < dxy, "dxy", muon_dxy);
    const double muon_dz = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muon_dz < dz, "dz", muon_dz);
    cut(muon->isTightMuon(*primaryVertex), "muonID");
    if(productionMode == ProductionMode::hh)
        cut(muon.GetIsolation() < pfRelIso04, "iso", muon.GetIsolation());
}

void TupleProducer_muTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::MuTau::tauID;

    cut(true, "gt0_tau_cand");
    const LorentzVector& p4 = tau.GetMomentum();
    const double pt_cut= productionMode == ProductionMode::h_tt_mssm ? cuts::H_tautau_2016_mssm::MuTau::tauID::pt : pt;
    cut(p4.Pt() > pt_cut, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < eta, "eta", p4.Eta());
    const auto dmFinding = tau->tauID("decayModeFinding");
    cut(dmFinding > decayModeFinding, "decayMode", dmFinding);
    auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
    if(productionMode == ProductionMode::hh) {
        cut(tau->tauID("againstElectronVLooseMVA6") > againstElectronVLooseMVA6, "againstElectron");
        cut(tau->tauID("againstMuonTight3") > againstMuonTight3, "againstMuon");
    }
}

void TupleProducer_muTau::FillEventTuple(const SelectionResults& selection)
{
    using namespace analysis;
    static constexpr float default_value = ntuple::DefaultFillValue<Float_t>();
    static constexpr int default_int_value = ntuple::DefaultFillValue<Int_t>();

    BaseTupleProducer::FillEventTuple(selection);
    eventTuple().channelID = static_cast<int>(analysis::Channel::MuTau);

    // Leg 1, lepton
    const MuonCandidate& muon = selection.higgs->GetFirstDaughter();
    eventTuple().p4_1     = LorentzVectorM(muon.GetMomentum());
    eventTuple().q_1      = muon.GetCharge();
    eventTuple().d0_1     = muon->muonBestTrack()->dxy(primaryVertex->position());
    eventTuple().dZ_1     = muon->muonBestTrack()->dz(primaryVertex->position());
    eventTuple().iso_1    = muon.GetIsolation();
    eventTuple().id_e_mva_nt_loose_1 = default_value;
    eventTuple().gen_match_1 = isMC ? gen_truth::genMatch(muon->p4(), *genParticles) : default_int_value;

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_muTau);
