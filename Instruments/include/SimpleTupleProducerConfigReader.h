/*! Definition of a reader for the simple tuple producer configuration.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SimpleTupleProducerConfig.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"

namespace simple_tuple {

class ProducerSetupEntryReader : public analysis::ConfigEntryReaderT<ProducerSetup> {
public:
    using Condition = ConfigEntryReader::Condition;
    using Base = ConfigEntryReaderT<ProducerSetup>;
    using Base::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("int_lumi", 1, Condition::equal_to);

        Base::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("int_lumi", current.int_lumi);
        ParseMappedEntryList("trigger_paths", current.trigger_paths, false);
    }
};

class SampleDescriptorEntryReader : public analysis::ConfigEntryReaderT<SampleDescriptor, SampleDescriptorCollection> {
public:
    using Condition = ConfigEntryReader::Condition;
    using Base = ConfigEntryReaderT<SampleDescriptor, SampleDescriptorCollection>;
    using Base::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("file_path", 1, Condition::equal_to);
        CheckReadParamCounts("channels", 1, Condition::less_equal);
        CheckReadParamCounts("sample_type", 1, Condition::equal_to);
        CheckReadParamCounts("higgs_prod", 1, Condition::less_equal);
        CheckReadParamCounts("higgs_gen_mass", 1, Condition::less_equal);
        CheckReadParamCounts("cross_section", 1, Condition::less_equal);

        Base::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("file_path", current.file_path);
        ParseEntryList("channels", current.channels);
        ParseEntry("sample_type", current.sample_type);
        ParseEntry("higgs_prod", current.higgs_prod);
        ParseEntry("higgs_gen_mass", current.higgs_gen_mass);
        ParseEntry<double, analysis::NumericalExpression>("cross_section", current.cross_section,
                                                          [](double xs){ return xs > 0; });
    }
};

} // namespace simple_tuple
