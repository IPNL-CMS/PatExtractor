#pragma once

namespace edm {
  class EventSetup;
  class Event;
}

class MCExtractor;

class SuperBaseExtractor
{
  public:
    virtual ~SuperBaseExtractor() {}
    virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor) = 0;
    virtual void getInfo(int ievt) = 0;

    bool isOK() {
      return m_OK;
    }

  protected:
    bool m_OK;

};

