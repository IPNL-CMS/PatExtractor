#pragma once

namespace edm {
  class EventSetup;
  class Event;
}

class MCExtractor;

class SuperBaseExtractor
{
  public:
    SuperBaseExtractor(): m_isMC(false), m_OK(false) {} 
    virtual ~SuperBaseExtractor() {}
    virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor) = 0;
    virtual void getInfo(int ievt) = 0;

    bool isOK() {
      return m_OK;
    }

    void setIsMC(bool isMC) {
      m_isMC = isMC;
    }

  protected:
    bool m_isMC;
    bool m_OK;
};

