#pragma once

/**
 * BaseExtractor
 * \brief: Base class for extracting info
 */


//Include PAT info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/View.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

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

template<class ObjectType>
class BaseExtractor: public SuperBaseExtractor
{

  public:
    BaseExtractor(const std::string& name):
      m_name(name) {}
    virtual ~BaseExtractor() {}

    virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor) {

      edm::Handle<edm::View<ObjectType>> handle;
      event.getByLabel(m_tag, handle);
      edm::View<ObjectType> collection = *handle;

      reset();
      fillSize(collection.size());

      int i = 0;
      for (auto& object: collection) {
        writeInfo(object, i);
        if (mcExtractor) {
          doMCMatch(object, mcExtractor, i);
        }
        i++;
      }

      fillTree();
    }

    virtual void writeInfo(const ObjectType& object, int index) = 0;
    virtual void doMCMatch(const ObjectType& object, MCExtractor* mcExtractor, int index) = 0;

    virtual void reset() = 0;
    virtual void fillTree() = 0;
    virtual void getInfo(int ievt) = 0;

    void fillSize(size_t size) {
      m_size = size;
    }

    int getSize() {
      return m_size;
    }

    void setPF(bool isPF) {
      m_isPF = isPF;
    }

  protected:

    std::string m_name;
    edm::InputTag m_tag;
    bool  m_isPF;
    TFile* m_file;

    size_t   m_size;
};
