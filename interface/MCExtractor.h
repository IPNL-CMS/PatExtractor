#pragma once

/**
 * MCExtractor
 * \brief: Base class for extracting MC info
 */


//Include PAT info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "../interface/SuperBaseExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class MCExtractor: public SuperBaseExtractor
{

 public:

  MCExtractor(const std::string& name, bool doTree, bool doJpsi = false);
  MCExtractor(const std::string& name, TFile *a_file, bool doJpsi = false);
  virtual ~MCExtractor();

  virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor);

  void reset();
  void fillTree(); 
  void getInfo(int ievt); 

  int getSize() const;

  // Setters/Getters

  int   getStatus(int index){return m_MC_status[index];}
  int   getType(int index){return m_MC_type[index];}
  float getPx(int index){return m_MC_px[index];}
  float getPy(int index){return m_MC_py[index];}
  float getPz(int index){return m_MC_pz[index];}
  float getE(int index){return m_MC_E[index];}
  float getMom1Index(int index){return m_MC_imot1[index];}
  TLorentzVector* getP4(int index) {
    return static_cast<TLorentzVector*>((*m_MC_lorentzvector)[index]);
  }

  int getPatIndex(int index) const { return m_MC_index[index]; }

 private:
  
  TTree* m_tree_MC;

  static const int 	m_MCs_MAX        = 1000;

  int   m_n_MCs;
  TClonesArray *m_MC_lorentzvector;
  int   m_MC_index[m_MCs_MAX];
  int   m_MC_status[m_MCs_MAX];
  int   m_MC_type[m_MCs_MAX];
  int   m_MC_imot1[m_MCs_MAX];
  int   m_MC_imot2[m_MCs_MAX];
  int   m_MC_generation[m_MCs_MAX];
  float m_MC_E[m_MCs_MAX];
  float	m_MC_px[m_MCs_MAX];
  float	m_MC_py[m_MCs_MAX];
  float	m_MC_pz[m_MCs_MAX];
  float	m_MC_vx[m_MCs_MAX];
  float	m_MC_vy[m_MCs_MAX];
  float	m_MC_vz[m_MCs_MAX];
  float	m_MC_eta[m_MCs_MAX];
  float	m_MC_phi[m_MCs_MAX];
  
  bool _doJpsi;
  bool m_MC_JPsiFromTop[m_MCs_MAX];
  bool m_MC_JPsiFromAntiTop[m_MCs_MAX];
  bool m_MC_LeptonFromTop[m_MCs_MAX];
  bool m_MC_LeptonFromAntiTop[m_MCs_MAX];

  void constructGeneration(int gene, int npart);
};
