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

    PhotonExtractor(const std::string& name, const edm::ParameterSet&);
    PhotonExtractor(const std::string& name, const edm::ParameterSet&, TFile *a_file);
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
    
    bool  passLooseId(int phidx) const { return m_pho_passLooseId[phidx]; }
    bool  passMediumId(int phidx) const { return m_pho_passMediumId[phidx]; }
    bool  passTightId(int phidx) const { return m_pho_passTightId[phidx]; }
    
    // Setters/Getters

  private:

    TTree* m_tree_photon;
    edm::InputTag m_photonTag;
    edm::InputTag m_rhoTag;
    edm::InputTag m_phoLooseIdMapTag;
    edm::InputTag m_phoMediumIdMapTag;
    edm::InputTag m_phoTightIdMapTag;

    edm::EDGetTokenT<pat::PhotonCollection> m_token;
    edm::EDGetTokenT<double> m_rhoToken;
    // ID decision objects
    edm::EDGetTokenT<edm::ValueMap<bool> > m_phoLooseIdMapToken;
    edm::EDGetTokenT<edm::ValueMap<bool> > m_phoMediumIdMapToken;
    edm::EDGetTokenT<edm::ValueMap<bool> > m_phoTightIdMapToken;

    static const int 	m_photons_MAX    = 100;

    float m_deltaR_cut;

    TClonesArray* m_pho_lorentzvector;
    float       m_pho_vx[m_photons_MAX];
    float	m_pho_vy[m_photons_MAX];
    float	m_pho_vz[m_photons_MAX];
    bool	m_pho_passLooseId[m_photons_MAX];
    bool	m_pho_passMediumId[m_photons_MAX];
    bool	m_pho_passTightId[m_photons_MAX];
    float	m_pho_chargedHadronsIsolation[m_photons_MAX];
    float	m_pho_neutralHadronsIsolation[m_photons_MAX];
    float	m_pho_photonIsolation[m_photons_MAX];
    int   m_pho_MCIndex[m_photons_MAX];
};

#endif 
