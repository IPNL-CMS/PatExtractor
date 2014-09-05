#pragma once

#include <memory>
#include <utility>

#include <TFile.h>

#include <Extractors/PatExtractor/interface/ScaleFactorService.h>

#include <FWCore/Framework/interface/ConsumesCollector.h>
#include <FWCore/PluginManager/interface/PluginFactory.h>

namespace edm {
  class EventSetup;
  class Event;
}

class MCExtractor;

class SuperBaseExtractor
{
  public:
    SuperBaseExtractor(const std::string& name, const edm::ParameterSet&): m_name(name), m_isMC(false), m_OK(false) {} 
    SuperBaseExtractor(const std::string& name, const edm::ParameterSet&, TFile *file): m_name(name), m_isMC(false), m_OK(false) {}
    virtual ~SuperBaseExtractor() {}

    virtual void doConsumes(edm::ConsumesCollector&& collector) {
      m_superCalled = true;
    }

    virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor) = 0;
    virtual void getInfo(int ievt) = 0;

    bool isOK() {
      return m_OK;
    }

    void setIsMC(bool isMC) {
      m_isMC = isMC;
    }

    virtual void beginJob(bool isInAnalysisMode) {}

    virtual void endJob(bool isInAnalysisMode) {}

    void check() const {
      if (! m_superCalled) {
        throw new std::logic_error("Base method 'SuperBaseExtractor::beginJob' was not called. Please check the extractors.");
      }
    }

    std::string getName() const {
      return m_name;
    }

  protected:
    std::string m_name;
    bool m_isMC;
    bool m_OK;

  private:
    bool m_superCalled = false;
    
};

typedef edmplugin::PluginFactory<SuperBaseExtractor* (const std::string&, const edm::ParameterSet&)> PatExtractorExtractorFactory;
typedef edmplugin::PluginFactory<SuperBaseExtractor* (const std::string&, const edm::ParameterSet&, TFile*)> PatExtractorExtractorReadOnlyFactory;
