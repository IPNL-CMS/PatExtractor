#ifndef MUONEXTRACTOR_H
#define MUONEXTRACTOR_H

/**
 * MuonExtractor
 * \brief: Base class for extracting muon info
 */


//Include PAT info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/View.h"

#include "../interface/BaseExtractor.h"
#include "../interface/MCExtractor.h"
#include <Extractors/PatExtractor/interface/ScaleFactor.h>

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class MuonExtractor: public BaseExtractor<pat::Muon>
{

  public:

    MuonExtractor(const std::string& name, const edm::InputTag& tag, const edm::InputTag& vertexTag, bool doTree);
    MuonExtractor(const std::string& name, TFile *a_file);
    virtual ~MuonExtractor();

    virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Muon& object, int index);

    void reset();
    void fillTree(); 
    void getInfo(int ievt); 

    virtual const reco::Candidate* getGenParticle(const pat::Muon& muon) {
      return muon.genLepton();
    }

    virtual void setGenParticleIndex(int genParticleIndex, int index) {
      m_muo_MCIndex[index] = genParticleIndex;
    }
    
    virtual float getMCMatchDeltaR() {
      return 0.5;
    }

    virtual float getMCMatchDPtRel() {
      return 0.5;
    }

    virtual std::vector<int> getPdgIds() {
      return {13};
    }

    virtual TLorentzVector getP4(const pat::Muon& object) {
      TLorentzVector p4;
      p4.SetPxPyPzE(object.px(), object.py(), object.pz(), object.energy());

      return p4;
    }

    // Setters/Getters

    TLorentzVector *getMuLorentzVector(int muidx) {return (TLorentzVector*)m_muo_lorentzvector->At(muidx);}
    int getMuisGlobal(int muidx) {return m_muo_isGlobal[muidx];}
    int getMuisTracker(int muidx) {return m_muo_isTracker[muidx];}
    int getMunValPixelHits(int muidx) {return m_muo_nValPixelHits[muidx];}
    int getMunValTrackerHits(int muidx) {return m_muo_nValTrackerHits[muidx];}
    int getMunMatches(int muidx) {return m_muo_nMatches[muidx];}
    float getMunormChi2(int muidx) {return m_muo_normChi2[muidx];}
    float getMudB(int muidx) {return m_muo_dB[muidx];}
    float getMupfChargedHadronIso(int muidx) {return m_muo_pfChargedHadronIso[muidx];}
    float getMupfNeutralHadronIso(int muidx) {return m_muo_pfNeutralHadronIso[muidx];}
    float getMupfPhotonIso(int muidx) {return m_muo_pfPhotonIso[muidx];}
    int getMuMCIndex(int muidx){return m_muo_MCIndex[muidx];}
    int getMuCharge(int muidx){return m_muo_charge[muidx];}

    float getTrackerLayersWithMeasurements(int index) {
      return m_muo_trackerLayersWithMeasurement[index];
    }

    float getGlobalTrackNumberOfValidMuonHits(int index) {
      return m_muo_globalTrackNumberOfValidHits[index];
    }

    int getNumberOfMatchedStations(int index) const {
      return m_muo_nMatchedStations[index];
    }

    float getdZ(int index) const {
      return m_muo_dZ[index];
    }

    float getDeltaBetaCorrectedRelativeIsolation(int index) const {
      return m_muo_deltaBetaCorrectedRelIsolation[index];
    }

    ScaleFactor getScaleFactor(int index) const {
      return m_scaleFactors.at(index);
    }


  private:

    edm::InputTag m_vertexTag;

    TTree* m_tree_muon;

    static const int 	m_muons_MAX  = 100;

    float m_deltaR_cut;

    TClonesArray* m_muo_lorentzvector;
    float	m_muo_vx[m_muons_MAX];
    float	m_muo_vy[m_muons_MAX];
    float	m_muo_vz[m_muons_MAX];
    int 	m_muo_isGlobal[m_muons_MAX];
    int 	m_muo_isTracker[m_muons_MAX];
    float m_muo_dB[m_muons_MAX];
    float m_muo_normChi2[m_muons_MAX];
    int 	m_muo_nValTrackerHits[m_muons_MAX]; // includes double sided layers and overlapped layers
    int 	m_muo_nValPixelHits[m_muons_MAX]; // without overlap
    int 	m_muo_nMatches[m_muons_MAX]; 
    float m_muo_trackIso[m_muons_MAX];
    float m_muo_ecalIso[m_muons_MAX];
    float m_muo_hcalIso[m_muons_MAX];
    int	m_muo_charge[m_muons_MAX];
    float m_muo_d0[m_muons_MAX];
    float m_muo_d0error[m_muons_MAX];
    float m_muo_pfParticleIso[m_muons_MAX]; // isolation calculated with all the PFCandidates

    float m_muo_pfChargedHadronIso[m_muons_MAX]; // isolation calculated with only the charged hadron PFCandidates
    float m_muo_pfNeutralHadronIso[m_muons_MAX]; // isolation calculated with only the neutral hadron PFCandidates
    float m_muo_pfPhotonIso[m_muons_MAX]; // Returns the isolation calculated with only the photon PFCandidates
    float m_muo_relIsolation[m_muons_MAX];
    float m_muo_deltaBetaCorrectedRelIsolation[m_muons_MAX];

    int   m_muo_MCIndex[m_muons_MAX];

    int   m_muo_nMatchedStations[m_muons_MAX];
    float m_muo_trackerLayersWithMeasurement[m_muons_MAX];
    float m_muo_dZ[m_muons_MAX];
    float m_muo_pixelLayerWithMeasurement[m_muons_MAX];
    float m_muo_globalTrackNumberOfValidHits[m_muons_MAX];

    ScaleFactorCollection m_scaleFactors;
};

#endif 
