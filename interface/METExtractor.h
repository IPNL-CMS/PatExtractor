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

#include "Extractors/PatExtractor/interface/BaseExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class METExtractor: public BaseExtractor<pat::MET>
{

 public:

  METExtractor(const std::string& name, const edm::InputTag& tag, bool doTree);
  METExtractor(const std::string& name, TFile *a_file);
  virtual ~METExtractor();


  virtual void writeInfo(const pat::MET& part, int index); 
  virtual void doMCMatch(const pat::MET& object, MCExtractor* mcExtractor, int index) {}

  void getInfo(int ievt);

  void reset();
  void fillTree(); 

  TLorentzVector *getMETLorentzVector(int metidx) {return (TLorentzVector*)m_met_lorentzvector->At(metidx);}

  void setMETLorentzVector(int idx, float E, float Px, float Py, float Pz);

 private:
  
  TTree* m_tree_met;
  static const int 	m_mets_MAX  = 2;
  TClonesArray* m_met_lorentzvector;
};

#endif 
