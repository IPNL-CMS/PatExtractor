#ifndef HLTEXTRACTOR_H
#define HLTEXTRACTOR_H

/**
 * HLTExtractor
 * \brief: Base class for extracting HLT info
 */

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

//Include std C++
#include <iostream>
#include <vector>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class HLTExtractor
{

 public:

  HLTExtractor(bool doTree);
  HLTExtractor(TFile *a_file);
  ~HLTExtractor();


  void writeInfo(const edm::Event *event); 
  void getInfo(int ievt); 
  void reset();
  void fillTree(); 
  void fillSize(int size);
  int  getSize();
  void print();

  // Setters/Getters

  bool isOK() {return m_OK;}

  std::string paths(int i) {return m_HLT_vector->at(i);}

 private:
  
  TTree* m_tree_HLT;

  int   m_n_HLTs;
  std::vector< std::string > *m_HLT_vector;

  bool m_OK;
};

#endif 
