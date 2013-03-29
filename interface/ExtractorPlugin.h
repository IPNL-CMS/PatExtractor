#pragma once

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

class PatExtractor;

namespace patextractor {

  class Plugin {
    public:
      Plugin(const edm::ParameterSet& config) {};
      
      void setIsMC(bool isMC) {
        m_isMC = isMC;
      }

      virtual void analyze(const edm::Event&, const edm::EventSetup&, PatExtractor&) = 0;
      virtual void analyze(const edm::EventSetup&, PatExtractor&) {};

    protected:
      bool m_isMC;
  };

}

typedef edmplugin::PluginFactory<patextractor::Plugin* (const edm::ParameterSet&)> PatExtractorPluginFactory;
