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

#include "../interface/MCExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class MuonExtractor
{

 public:

  MuonExtractor(bool doTree,edm::InputTag tag);
  MuonExtractor(TFile *a_file);
  ~MuonExtractor();

  void writeInfo(const edm::Event *event,MCExtractor* m_MC, bool doMC); 

  void writeInfo(const pat::Muon *part, int index); 

  void reset();
  void fillTree(); 
  void fillSize(int size);
  int  getSize();
  void setPF(bool isPF);
  void getInfo(int ievt); 

  int getMatch(const pat::Muon *part, MCExtractor* m_MC);

  // Setters/Getters

  bool isOK() {return m_OK;}

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

 private:
  
  TTree* m_tree_muon;

  static const int 	m_muons_MAX  = 100;

  edm::InputTag m_tag;
  float m_deltaR_cut;
  
  bool  m_isPF_muon;

  int   m_n_muons;
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
  int   m_muo_MCIndex[m_muons_MAX];

  bool m_OK;
};

#endif 
