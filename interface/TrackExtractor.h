#ifndef TRACKEXTRACTOR_H
#define TRACKEXTRACTOR_H

/**
 * TrackExtractor
 * \brief: Base class for extracting track info
 */


//Include PAT info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "../interface/BaseExtractor.h"
#include "../interface/MCExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"

class TrackExtractor: public BaseExtractor<reco::Track>
{

  public:

    TrackExtractor(const std::string& name, const edm::InputTag& tag, bool doTree);
    TrackExtractor(const std::string& name, TFile* a_file);
    virtual ~TrackExtractor();

    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const reco::Track& part, int index);

    virtual const reco::Candidate* getGenParticle(const reco::Track& track) {
      return NULL;
    }

    virtual void setGenParticleIndex(int genParticleIndex, int index) {
    }

    virtual float getMCMatchDeltaR() {
      return 0.;
    }

    virtual float getMCMatchDPtRel() {
      return 0.;
    }

    virtual std::vector<int> getPdgIds() {
      return {};
    }

    virtual TLorentzVector getP4(const reco::Track& object) {
      TLorentzVector p4;
      return p4;
    }

    void reset();
    void fillTree();
    void getInfo(int ievt); 

    // Setters/Getters

    float getTrackpx(int muidx) {return m_trk_px[muidx];}
    float getTrackpy(int muidx) {return m_trk_py[muidx];}
    float getTrackpz(int muidx) {return m_trk_pz[muidx];}
    float getMunormChi2(int muidx) {return m_trk_normChi2[muidx];}

  private:

    TTree* m_tree_track;

    static const int 	m_tracks_MAX  = 2000;

    float	m_trk_px[m_tracks_MAX];
    float	m_trk_py[m_tracks_MAX];
    float	m_trk_pz[m_tracks_MAX];
    float	m_trk_vx[m_tracks_MAX];
    float	m_trk_vy[m_tracks_MAX];
    float	m_trk_vz[m_tracks_MAX];
    float	m_trk_eta[m_tracks_MAX];
    float	m_trk_phi[m_tracks_MAX];
    float m_trk_normChi2[m_tracks_MAX];
    int 	m_trk_nValidHits[m_tracks_MAX]; 
    int	m_trk_charge[m_tracks_MAX];
    float m_trk_d0[m_tracks_MAX];

};

#endif 
