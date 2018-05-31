/*! Definition of a simple tuple producer configuration.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <unordered_map>
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/map_vec.h"
#include "h-tautau/Analysis/include/SimpleTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"

namespace simple_tuple {

struct ProducerSetup {
    std::string name;
    double int_lumi{0};
    std::map<analysis::Channel, std::vector<std::string>> trigger_paths;
};

using ProducerSetupCollection = std::unordered_map<std::string, ProducerSetup>;

struct SampleDescriptor {
    std::string name;
    std::string file_path;
    std::set<analysis::Channel> channels;
    SampleType sample_type{SampleType::Data};
    HiggsProductionType higgs_prod{HiggsProductionType::Undefined};
    int higgs_gen_mass{0};
    double cross_section{1};
};

using SampleDescriptorCollection = analysis::map_vec<std::string, SampleDescriptor>;

} // namespace simple_tuple

