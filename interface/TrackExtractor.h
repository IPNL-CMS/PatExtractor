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



//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"

class TrackExtractor
{

 public:

  TrackExtractor(bool doTree,edm::InputTag tag);
  TrackExtractor(TFile* a_file);
  ~TrackExtractor();

  void writeInfo(const edm::Event *event); 

  void writeInfo(const reco::Track *part, int index); 

  void reset();
  void fillTree(); 
  void fillSize(int size);
  int  getSize();
  void getInfo(int ievt); 

  // Setters/Getters

  bool isOK() {return m_OK;}
  float getTrackpx(int muidx) {return m_trk_px[muidx];}
  float getTrackpy(int muidx) {return m_trk_py[muidx];}
  float getTrackpz(int muidx) {return m_trk_pz[muidx];}
  float getMunormChi2(int muidx) {return m_trk_normChi2[muidx];}

 private:
  
  TTree* m_tree_track;

  static const int 	m_tracks_MAX  = 2000;

  edm::InputTag m_tag;


  bool m_OK;

  int   m_n_tracks;
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
