#ifndef PHOTONEXTRACTOR_H
#define PHOTONEXTRACTOR_H

/**
 * PhotonExtractor
 * \brief: Base class for extracting photon info
 */


//Include PAT info
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "../interface/BaseExtractor.h"
#include "../interface/MCExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class PhotonExtractor: public BaseExtractor<pat::Photon>
{

  public:

    PhotonExtractor(const std::string& name, const edm::InputTag& tag, bool doTree);
    PhotonExtractor(const std::string& name, TFile *a_file);
    virtual ~PhotonExtractor();
    virtual void doConsumes(edm::ConsumesCollector&& collector);
    
    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* m_MC); 

    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Photon& part, int index, const pat::PhotonRef& photonRef);
    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Photon& part, int index) {} 
    
    void getInfo(int ievt);

    void reset();
    void fillTree(); 

    virtual const reco::Candidate* getGenParticle(const pat::Photon& photon) {
      return photon.genPhoton();
    }

    virtual void setGenParticleIndex(int genParticleIndex, int index) {
      m_pho_MCIndex[index] = genParticleIndex;
    }

    virtual float getMCMatchDeltaR() {
      return 0.2;
    }

    virtual float getMCMatchDPtRel() {
      return 1.0;
    }

    virtual std::vector<int> getPdgIds() {
      return {22};
    }

    virtual TLorentzVector getP4(const pat::Photon& object) {
      TLorentzVector p4;
      p4.SetPxPyPzE(object.px(), object.py(), object.pz(), object.energy());

      return p4;
    }
    
    TLorentzVector *getP4(int photidx) {return (TLorentzVector*)m_pho_lorentzvector->At(photidx);}
    
    bool  hasPixelSeed(int phidx) const { return m_pho_hasPixelSeed[phidx]; }
    
    float getHadTowOverEm(int phidx) const { return m_pho_hadTowOverEm[phidx]; }
    float getSigmaIetaIeta(int phidx) const { return m_pho_sigmaIetaIeta[phidx]; }
    bool  hasMatchedPromptElectron(int phidx) const { return m_pho_hasMatchedPromptElectron[phidx]; }
    float getChargedHadronsIsolation(int phidx) const { return m_pho_chargedHadronsIsolation[phidx]; }
    float getNeutralHadronsIsolation(int phidx) const { return m_pho_neutralHadronsIsolation[phidx]; }
    float getPhotonIsolation(int phidx) const { return m_pho_photonIsolation[phidx]; }
    
    // Setters/Getters

  private:

    TTree* m_tree_photon;
    edm::EDGetTokenT<double> m_rhoToken;
    edm::EDGetTokenT<edm::ValueMap<bool>> m_matchedPromptElectronToken;
    edm::EDGetTokenT<edm::ValueMap<double>> m_chargedHadronsIsolationToken;
    edm::EDGetTokenT<edm::ValueMap<double>> m_neutralHadronsIsolationToken;
    edm::EDGetTokenT<edm::ValueMap<double>> m_photonsIsolationToken;

    static const int 	m_photons_MAX    = 100;

    float m_deltaR_cut;

    TClonesArray* m_pho_lorentzvector;
    float       m_pho_vx[m_photons_MAX];
    float	m_pho_vy[m_photons_MAX];
    float	m_pho_vz[m_photons_MAX];
    bool	m_pho_hasPixelSeed[m_photons_MAX];
    float	m_pho_hadTowOverEm[m_photons_MAX];
    float	m_pho_sigmaIetaIeta[m_photons_MAX];
    bool	m_pho_hasMatchedPromptElectron[m_photons_MAX];
    float	m_pho_chargedHadronsIsolation[m_photons_MAX];
    float	m_pho_neutralHadronsIsolation[m_photons_MAX];
    float	m_pho_photonIsolation[m_photons_MAX];
    int   m_pho_MCIndex[m_photons_MAX];
};

#endif 
