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
    
    /**
     * \brief Declares which data an inheriting plugin is going to read from an EDM event
     * 
     * See documentation in [1].
     * [1] https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/FWCore/Framework/interface/ConsumesCollector.h
     */
    virtual void doConsumes(edm::ConsumesCollector&& collector) {
      m_superCalled = true;
    }
    
    /// Processes current event in the source EDM file
    virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor) = 0;
    
    virtual void getInfo(int ievt) = 0;
    
    /// Reports if the extractor has been created without errors and can be used as intended
    bool isHealthy() const
    {
        return m_OK;
    }
    
    /**
     * \brief Proxy for isHealthy
     * 
     * \deprecated Use the isHealthy method instead
     */
    bool isOK() const __attribute__ ((deprecated))
    {
        return isHealthy();
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
    /// Specifies whether the extractor has been created without errors
    void setHealthy(bool ok)
    {
        m_OK = ok;
    }
    
protected:
    std::string m_name;
    bool m_isMC;
    
private:
    bool m_superCalled = false;
    
    /// Indicates if the extractor has been constructed without errors
    bool m_OK;
};

typedef edmplugin::PluginFactory<SuperBaseExtractor* (const std::string&, const edm::ParameterSet&)> PatExtractorExtractorFactory;
typedef edmplugin::PluginFactory<SuperBaseExtractor* (const std::string&, const edm::ParameterSet&, TFile*)> PatExtractorExtractorReadOnlyFactory;
