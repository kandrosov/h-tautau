/*! Definition of a tuple with minimal information needed to perform a simpified analysis.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"

#define SIMPLE_EVENT_DATA() \
    VAR(ULong64_t, event_id) /* Event identifier */ \
    VAR(Int_t, sample_type) /* Sample type */ \
    VAR(Int_t, higgs_prod_type) /* Higgs production mechanism */ \
    VAR(Int_t, higgs_gen_mass) /* Higgs mass at the generator level */ \
    VAR(Int_t, channel_id) /* Channel: eTau, muTau or tauTau */ \
    VAR(Float_t, weight) /* Event weight */ \
    VAR(Int_t, npv) /* Number of the primary vertices in the event */ \
    /* The first tau candidate (electron for eTau channel,  muon for muTau channel, pt leading hadronic tau in tauTau channel */ \
    VAR(Float_t, px_1) /* The x component of the 4-momentum of the first tau candidate */ \
    VAR(Float_t, py_1) /* The y component of the 4-momentum of the first tau candidate */ \
    VAR(Float_t, pz_1) /* The z component of the 4-momentum of the first tau candidate */ \
    VAR(Float_t, E_1) /* The energy component of the 4-momentum of the first tau candidate */ \
    VAR(Int_t, q_1) /* Charge of the first tau candidate */ \
    VAR(Float_t, dxy_1) /* dxy with respect to the primary vertex of the first tau candidate */ \
    VAR(Float_t, dz_1) /* dz with respect to the primary vertex of the first tau candidate */ \
    VAR(Float_t, decay_mode_1) /* Decay mode of the first tau candidate */ \
    VAR(Int_t, gen_match_1) /* Generator matching result for the first tau candidate */\
    VAR(Float_t, gen_px_1) /* The x component of the 4-momentum of the gen particle matched with the first tau candidate */ \
    VAR(Float_t, gen_py_1) /* The y component of the 4-momentum of the gen particle matched with the first tau candidate */ \
    VAR(Float_t, gen_pz_1) /* The z component of the 4-momentum of the gen particle matched with the first tau candidate */ \
    VAR(Float_t, gen_E_1) /* The energy component of the 4-momentum of the gen particle matched with the first tau candidate */ \
    /* The second tau candidate (hadronic tau for eTau and muTau channels, pt sub-leading hadronic tau in tauTau channel */ \
    VAR(Float_t, px_2) /* The x component of the 4-momentum of the second tau candidate */ \
    VAR(Float_t, py_2) /* The y component of the 4-momentum of the second tau candidate */ \
    VAR(Float_t, pz_2) /* The z component of the 4-momentum of the second tau candidate */ \
    VAR(Float_t, E_2) /* The energy component of the 4-momentum of the second tau candidate */ \
    VAR(Int_t, q_2) /* Charge of the second tau candidate */ \
    VAR(Float_t, dxy_2) /* dxy with respect to the primary vertex of the second tau candidate  */ \
    VAR(Float_t, dz_2) /* dz with respect to the primary vertex of the second tau candidate */ \
    VAR(Float_t, decay_mode_2) /* Decay mode of the second tau candidate */ \
    VAR(Int_t, gen_match_2) /* Generator matching result for the second tau candidate */\
    VAR(Float_t, gen_px_2) /* The x component of the 4-momentum of the gen particle matched with the second tau candidate */ \
    VAR(Float_t, gen_py_2) /* The y component of the 4-momentum of the gen particle matched with the second tau candidate */ \
    VAR(Float_t, gen_pz_2) /* The z component of the 4-momentum of the gen particle matched with the second tau candidate */ \
    VAR(Float_t, gen_E_2) /* The energy component of the 4-momentum of the gen particle matched with the second tau candidate */ \
    /* Missing transvers energy (MET) and the fitted mass of the Higgs candidate */ \
    VAR(Float_t, met_px) /* The x component of the MET */ \
    VAR(Float_t, met_py) /* The y component of the MET */ \
    VAR(Float_t, met_cov_00) /* (0, 0) element of the MET covariance matrix */ \
    VAR(Float_t, met_cov_01) /* (0, 1) element of the MET covariance matrix */ \
    VAR(Float_t, met_cov_10) /* (1, 0) element of the MET covariance matrix  */ \
    VAR(Float_t, met_cov_11) /* (1, 1) element of the MET covariance matrix  */ \
    VAR(Float_t, fitted_H_mass) /* Mass of the Higgs candidate obtained by the MCMC integration */ \
    /* Selected leading and sub-leading jet candidates in the event (excluding hadronic taus) */ \
    VAR(Float_t, n_jets) /* Total number of the selected jets in the event */ \
    VAR(Float_t, leading_jet_px) /* The x component of the 4-momentum of the leading jet */ \
    VAR(Float_t, leading_jet_py) /* The y component of the 4-momentum of the leading jet */ \
    VAR(Float_t, leading_jet_pz) /* The z component of the 4-momentum of the leading jet */ \
    VAR(Float_t, leading_jet_E) /* The energy component of the 4-momentum of the leading jet */ \
    VAR(Float_t, leading_jet_btag) /* Value of the b tagging discriminator for the leading jet */ \
    VAR(Float_t, subleading_jet_px) /* The x component of the 4-momentum of the sub-leading jet */ \
    VAR(Float_t, subleading_jet_py) /* The y component of the 4-momentum of the sub-leading jet */ \
    VAR(Float_t, subleading_jet_pz) /* The z component of the 4-momentum of the sub-leading jet */ \
    VAR(Float_t, subleading_jet_E) /* The energy component of the 4-momentum of the sub-leading jet */ \
    VAR(Float_t, subleading_jet_btag) /* Value of the b tagging discriminator for the sub-leading jet */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(simple_tuple, Event, EventTuple, SIMPLE_EVENT_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(simple_tuple, EventTuple, SIMPLE_EVENT_DATA)
#undef VAR
#undef SIMPLE_EVENT_DATA

namespace simple_tuple {
template<typename T>
constexpr T DefaultFillValue() { return T(-9e9); }

using ::analysis::operator<<;
using ::analysis::operator>>;

enum class SampleType { Data = 0, H_tautau = 1, DY = 2, TT = 3, W = 4, DiBoson = 5, SingleTop = 6, EWK = 7, QCD = 8 };
ENUM_NAMES(SampleType) = {
    { SampleType::Data, "Data" },
    { SampleType::H_tautau, "H_tautau" },
    { SampleType::DY, "DY" },
    { SampleType::TT, "TT" },
    { SampleType::W, "W" },
    { SampleType::QCD, "QCD" },
    { SampleType::DiBoson, "DiBoson" },
    { SampleType::SingleTop, "SingleTop" },
    { SampleType::EWK, "EWK" },
};

enum class HiggsProductionType { Undefined = 0, ggH = 1, VBF = 2, WH = 3, ZH = 4, ttH = 5 };
ENUM_NAMES(HiggsProductionType) = {
    { HiggsProductionType::Undefined, "Undefined" },
    { HiggsProductionType::ggH, "ggH" },
    { HiggsProductionType::VBF, "VBF" },
    { HiggsProductionType::WH, "WH" },
    { HiggsProductionType::ZH, "ZH" },
    { HiggsProductionType::ttH, "ttH" }
};

} // namespace simple_tuple
