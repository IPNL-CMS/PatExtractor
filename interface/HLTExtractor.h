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

#include "Extractors/PatExtractor/interface/BaseExtractor.h"
#include "Extractors/PatExtractor/interface/MCExtractor.h"

//Include std C++
#include <iostream>
#include <vector>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class HLTExtractor: public SuperBaseExtractor
{

 public:

  HLTExtractor(const std::string& name, bool doTree);
  HLTExtractor(const std::string& name, TFile *a_file);
  ~HLTExtractor();


  virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor);
  void getInfo(int ievt); 

  void reset();
  void fillTree(); 

  int getSize() const;

  // Setters/Getters

  std::string paths(int i) {return m_HLT_vector->at(i);}

  std::vector<std::string>* getPaths() {
    return m_HLT_vector;
  }

 private:
  
  TTree* m_tree_HLT;

  int m_n_HLTs;
  std::vector< std::string > *m_HLT_vector;

  bool m_OK;
};

#endif 
