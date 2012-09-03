#ifndef METEXTRACTOR_H
#define METEXTRACTOR_H

/**
 * METExtractor
 * \brief: Base class for extracting MET info
 */


//Include PAT info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/View.h"



//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class METExtractor
{

 public:

  METExtractor(bool doTree,edm::InputTag tag);
  METExtractor(TFile *a_file);
  ~METExtractor();


  void writeInfo(const pat::MET *part, int index); 
  void writeInfo(const edm::Event *event); 
  void getInfo(int ievt); 

  void reset();
  void fillTree(); 
  void fillSize(int size);
  int  getSize();
  TLorentzVector *getMETLorentzVector(int metidx) {return (TLorentzVector*)m_met_lorentzvector->At(metidx);}

  void setMETLorentzVector(int idx, float E, float Px, float Py, float Pz);

  bool isOK() {return m_OK;}
 private:
  
  TTree* m_tree_met;

  static const int 	m_mets_MAX  = 2;

  edm::InputTag m_tag;

  bool m_OK;
  int   m_n_mets;

  TClonesArray* m_met_lorentzvector;

};

#endif 
