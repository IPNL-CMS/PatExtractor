#ifndef ELECTRONEXTRACTOR_H
#define ELECTRONEXTRACTOR_H

/**
 * ElectronExtractor
 * \brief: Base class for extracting electron info
 */


//Include PAT info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/View.h"

#include "../interface/MCExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class ElectronExtractor
{

 public:

  ElectronExtractor(bool doTree,edm::InputTag tag);
  ElectronExtractor(TFile *a_file);
  ~ElectronExtractor();

  void writeInfo(const edm::Event *event,MCExtractor* m_MC, bool doMC); 
  void writeInfo(const pat::Electron *part, int index); 

  void reset();
  void fillTree(); 
  void fillSize(int size);
  int  getSize();
  void setPF(bool isPF); 
  void getInfo(int ievt); 
  
  int getMatch(const pat::Electron *part, MCExtractor* m_MC);

  bool isOK() {return m_OK;}
  TLorentzVector *getEleLorentzVector(int eidx) {return (TLorentzVector*)m_ele_lorentzvector->At(eidx);}
  float getEledB(int eidx) {return m_ele_dB[eidx];}
  float getElepfChargedHadronIso(int eidx) {return m_ele_pfChargedHadronIso[eidx];}
  float getElepfNeutralHadronIso(int eidx) {return m_ele_pfNeutralHadronIso[eidx];}
  float getElepfPhotonIso(int eidx) {return m_ele_pfPhotonIso[eidx];}

  int getEleMCIndex(int eidx){return m_ele_MCIndex[eidx];}
  int getEleCharge(int eidx) {return m_ele_charge[eidx];}
  int getEleHyperTight1MC(int eidx){return m_ele_eidHyperTight1MC[eidx];}
  int getEleLooseMC(int eidx)      {return m_ele_eidLooseMC[eidx];}
  int getEleMediumMC(int eidx)     {return m_ele_eidMediumMC[eidx];}
  int getEleSuperTightMC(int eidx) {return m_ele_eidSuperTightMC[eidx];}
  int getEleTightMC(int eidx)      {return m_ele_eidTightMC[eidx];}
  int getEleVeryLooseMC(int eidx)  {return m_ele_eidVeryLooseMC[eidx];}
  int getElepf_evspi(int eidx)     {return m_ele_eidpf_evspi[eidx];}
  int getElepf_evsmu(int eidx)     {return m_ele_eidpf_evsmu[eidx];}
  
 private:
  
  TTree* m_tree_electron;

  static const int 	m_electrons_MAX  = 100;


  edm::InputTag m_tag;
  float m_deltaR_cut;
  
  bool  m_isPF_electron;

  bool m_OK;
  int   m_n_electrons;

  TClonesArray* m_ele_lorentzvector;
  float	m_ele_vx[m_electrons_MAX];
  float	m_ele_vy[m_electrons_MAX];
  float	m_ele_vz[m_electrons_MAX];
  int	m_ele_charge[m_electrons_MAX];

  // electron id's
  
  /// old one
  int   m_ele_eidLoose[m_electrons_MAX]; 
  int   m_ele_eidRobustHighEnergy[m_electrons_MAX]; 
  int   m_ele_eidRobustLoose[m_electrons_MAX]; 
  int   m_ele_eidRobustTight[m_electrons_MAX]; 
  int   m_ele_eidTight[m_electrons_MAX]; 
  int   m_ele_eidpf_evspi[m_electrons_MAX]; 
  int   m_ele_eidpf_evsmu[m_electrons_MAX];
  
  /// new one
  int m_ele_eidHyperTight1MC[m_electrons_MAX];
  int m_ele_eidLooseMC[m_electrons_MAX];
  int m_ele_eidMediumMC[m_electrons_MAX];
  int m_ele_eidSuperTightMC[m_electrons_MAX];
  int m_ele_eidTightMC[m_electrons_MAX];
  int m_ele_eidVeryLooseMC[m_electrons_MAX];

  
  //
  int eidBit;
  bool pass ;
  float m_ele_dB[m_electrons_MAX];
  float m_ele_trackIso[m_electrons_MAX];
  float m_ele_ecalIso[m_electrons_MAX];
  float m_ele_hcalIso[m_electrons_MAX];
  float m_ele_pfParticleIso[m_electrons_MAX]; // isolation calculated with all the PFCandidates
  float m_ele_pfChargedHadronIso[m_electrons_MAX]; // isolation calculated with only the charged hadron PFCandidates
  float m_ele_pfNeutralHadronIso[m_electrons_MAX]; // isolation calculated with only the neutral hadron PFCandidates
  float m_ele_pfPhotonIso[m_electrons_MAX]; // Returns the isolation calculated with only the photon PFCandidates
  int   m_ele_numberOfMissedInnerLayer[m_electrons_MAX]; // Access the hit pattern counting (in the Tracker) the number of expected crossed layers  before the first trajectory's hit
  int   m_ele_MCIndex[m_electrons_MAX];
  

};

#endif 
