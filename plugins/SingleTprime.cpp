#include "SingleTprime.h"

SingleTprime::SingleTprime(const edm::ParameterSet& iConfig): Plugin(iConfig)
{
  // Initialize the analysis parameters using the ParameterSet iConfig
  int an_option = iConfig.getUntrackedParameter<int>("an_option", 0);
}

SingleTprime::::analysis(const edm::EventSetup& iSetup, PatExtractor& extractor)
{
  // Do the analysis
}

// Register the plugin inside the factory
DEFINE_EDM_PLUGIN(PatExtractorPluginFactory, SingleTprime,  "SingleTprime");
