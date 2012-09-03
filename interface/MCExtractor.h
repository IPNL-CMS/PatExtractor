#ifndef MCEXTRACTOR_H
#define MCEXTRACTOR_H

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


//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class MCExtractor
{

 public:

  MCExtractor(bool doTree);
  MCExtractor(TFile *a_file);
  ~MCExtractor();


  void writeInfo(const edm::Event *event); 

  void reset();
  void fillTree(); 
  void fillSize(int size);
  int  getSize();
  void getInfo(int ievt); 

  // Setters/Getters

  bool isOK() {return m_OK;}

  int   getStatus(int index){return m_MC_status[index];}
  int   getType(int index){return m_MC_type[index];}
  float getPx(int index){return m_MC_px[index];}
  float getPy(int index){return m_MC_py[index];}
  float getPz(int index){return m_MC_pz[index];}
  float getE(int index){return m_MC_E[index];}
  float getMom1Index(int index){return m_MC_imot1[index];}


 private:
  
  TTree* m_tree_MC;

  static const int 	m_MCs_MAX        = 1000;

  bool m_OK;

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

  void constructGeneration(int gene, int npart);
};

#endif 
