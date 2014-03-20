#pragma once

#include <memory>
#include <Extractors/PatExtractor/interface/ScaleFactorService.h>

namespace edm {
  class EventSetup;
  class Event;
}

class MCExtractor;

class SuperBaseExtractor
{
  public:
    SuperBaseExtractor(std::shared_ptr<ScaleFactorService> sf = std::shared_ptr<ScaleFactorService>()): m_isMC(false), m_OK(false), m_scaleFactorService(sf) {} 
    virtual ~SuperBaseExtractor() {}
    virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor) = 0;
    virtual void getInfo(int ievt) = 0;

    bool isOK() {
      return m_OK;
    }

    void setIsMC(bool isMC) {
      m_isMC = isMC;
    }

    virtual void beginJob() {}
    virtual void endJob() {}

  protected:
    bool m_isMC;
    bool m_OK;

    std::shared_ptr<ScaleFactorService> m_scaleFactorService;
    
};

