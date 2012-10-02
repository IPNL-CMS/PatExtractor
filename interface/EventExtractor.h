#ifndef EVENTEXTRACTOR_H
#define EVENTEXTRACTOR_H

/**
 * EventExtractor
 * \brief: Base class for extracting Event info
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "../interface/BaseExtractor.h"
#include "../interface/MCExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class EventExtractor: public SuperBaseExtractor
{

 public:

  EventExtractor(const std::string& name);
  EventExtractor(const std::string& name, TFile *a_file);
  virtual ~EventExtractor();


  virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor);

  void getInfo(int ievt);

  void reset();
  void print();

  // Setters/Getters

  int evtID() {return m_evtID;}
  int BCID() {return m_BCID;}
  int time() {return m_time;}
  int lumi() {return m_lumi;}
  int run() {return m_run;}
  int n_events() {return m_n_events;}
  int nPU() {return m_nPU;}
  float nTrueInteractions() { return m_nTrueInteractions; }
 private:
  
  TTree* m_tree_event;

  unsigned int   m_evtID;
  int   m_BCID;
  unsigned int   m_time;
  unsigned int   m_lumi;
  unsigned int   m_run; 
  unsigned int   m_n_events;
  int   m_nPU;
  float m_nTrueInteractions;
};

#endif 
