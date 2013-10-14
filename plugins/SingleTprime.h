#include <Extractors/PatExtractor/interface/ExtractorPlugin.h>

class SingleTprime: patextractor::Plugin {
  public:
    SingleTprime(const edm::ParameterSet& iConfig);

    virtual void analyze(const edm::EventSetup& iSetup, PatExtractor& extractor);
};
