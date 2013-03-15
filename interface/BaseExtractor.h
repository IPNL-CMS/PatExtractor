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

#include "DataFormats/Candidate/interface/Candidate.h"
#include "Extractors/PatExtractor/interface/MCExtractor.h"
#include "Extractors/PatExtractor/interface/SuperBaseExtractor.h"

//Include std C++
#include <iostream>
#include <limits>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

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
      if (!handle.isValid())
        return;

      edm::View<ObjectType> collection = *handle;

      reset();
      fillSize(collection.size());

      int i = 0;
      for (auto& object: collection) {
        writeInfo(event, iSetup, object, i);
        if (mcExtractor) {
          doMCMatch(object, event, mcExtractor, i);
        }
        i++;
      }

      fillTree();
    }

    virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const ObjectType& object, int index) = 0;

    virtual const reco::Candidate* getGenParticle(const ObjectType& object) = 0;
    virtual void setGenParticleIndex(int genParticleIndex, int index) = 0;

    virtual float getMCMatchDeltaR() = 0;
    virtual float getMCMatchDPtRel() = 0;
    virtual std::vector<int> getPdgIds() = 0;
    virtual TLorentzVector getP4(const ObjectType& object) = 0;

    virtual void reset() = 0;
    virtual void fillTree() = 0;
    virtual void getInfo(int ievt) = 0;

    void fillSize(uint32_t size) {
      m_size = size;
    }

    uint32_t getSize() {
      return m_size;
    }

    void setPF(bool isPF) {
      m_isPF = isPF;
    }

  protected:
    virtual void doMCMatch(const ObjectType& object, const edm::Event& event, MCExtractor* mcExtractor, int index) {

      // For now, always redo mc matching
      if (getPdgIds().size() == 0) {
        setGenParticleIndex(-4, index);
        return;
      }

      findMCMatch(object, event, mcExtractor, index);
      return;

      const reco::Candidate* gen = getGenParticle(object);
      if (gen == NULL) {
        if (getPdgIds().size() == 0) {
          setGenParticleIndex(-4, index);
        } else {
          findMCMatch(object, event, mcExtractor, index);
        }
        return;
      }

      if (false) {
        std::cout << "Matching MC for " << gen->pdgId() << std::endl;
        std::cout << "P: " << gen->px() << " " << gen->py() << " " << gen->pz() << std::endl;
      }

      if (gen->status() != 3) {
        findMCMatch(object, event, mcExtractor, index);
        return;
      }

      for (int i = 0; i < mcExtractor->getSize(); i++) {

        if (gen->pdgId() != mcExtractor->getType(i))
          continue;

        if (fabs(gen->px() - mcExtractor->getPx(i)) > 0.0001)
          continue;

        if (fabs(gen->py() - mcExtractor->getPy(i)) > 0.0001)
          continue;

        if (fabs(gen->pz() - mcExtractor->getPz(i)) > 0.0001)
          continue;

        setGenParticleIndex(i, index);
        return;
      }

      setGenParticleIndex(-6, index);
    }

    void findMCMatch(const ObjectType& object, const edm::Event& event, MCExtractor* mcExtractor, int index) {

      float dR_min = std::numeric_limits<float>::infinity();
      int mcIndex = -10;

      float dR_limit = getMCMatchDeltaR();
      float dP_limit = getMCMatchDPtRel();
      std::vector<int> pdgIds = getPdgIds();

      TLorentzVector p4 = getP4(object);

      int size = mcExtractor->getSize();
      for (int i = 0; i < size; i++) {
        // Check pdg id
        int pdgId = (int) fabs(mcExtractor->getType(i));
        if (std::find(pdgIds.begin(), pdgIds.end(), pdgId) == pdgIds.end()) {
          // Id not found, continue
          continue;
        }

        // Compute dR and dP, and check limits
        TLorentzVector mcp4; // = mcExtractor->getP4(i);
        mcp4.SetPxPyPzE(mcExtractor->getPx(i), mcExtractor->getPy(i), mcExtractor->getPz(i), mcExtractor->getE(i));
        float dR = mcp4.DeltaR(p4);
        if (dR > dR_limit)
          continue;

        float dP = fabs(mcp4.Pt() - p4.Pt()) / p4.Pt();
        if (dP > dP_limit)
          continue;

        if (dR < dR_min) {
          dR_min = dR;
          mcIndex = i;
        }
      }

      setGenParticleIndex(mcIndex, index);
      if (false) {
        std::cout << "\tType: " << mcExtractor->getType(mcIndex) << std::endl;
        std::cout << "\tP: " << mcExtractor->getPx(mcIndex) << " " << mcExtractor->getPy(mcIndex) << " " << mcExtractor->getPz(mcIndex) << std::endl;
      }
    }

    std::string m_name;
    edm::InputTag m_tag;

    bool  m_isPF;
    TFile* m_file;

    uint32_t m_size;
};
