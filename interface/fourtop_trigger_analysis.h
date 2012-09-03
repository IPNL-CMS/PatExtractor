#ifndef FOURTOP_TRIGGER_ANALYSIS_H
#define FOURTOP_TRIGGER_ANALYSIS_H

/*****************************

Simple example class showing how to perform an 
analysis using the PatExtractor tools

S.Viret (viret@in2p3.fr): 31/05/11

More info: http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.PHYTuto

 ******************************/

//Include std C++
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
using namespace std;

#include "AnalysisSettings.h"
#include "HLTExtractor.h"
#include "MCExtractor.h"
#include "TLorentzVector.h"

class fourtop_trigger_analysis
{
 public:
  fourtop_trigger_analysis(AnalysisSettings *settings);

  ~fourtop_trigger_analysis();
  
  //Selection

  //  bool isMuSel(TLorentzVector *muon);
  int  fourtop_trigger_Sel(HLTExtractor *hlt, MCExtractor *mc, int evtnum) ;

  void fourtop_trigger_finalize(int evtnum);
  void print_results(int index,std::string pathname,int evtnum);
  int  nSSlept(MCExtractor *mc);

  void reset();
  void fillTree();
   
 private:

  TTree* m_tree_fourtop_trigger;

  int n_tot_evt;
  int m_tot[500];             // Total number of events for a given HLT path
  std::multimap< std::string, int > m_HLT_vector_tot;  // Map of all the HLT paths with the number of times it fired
  std::multimap< std::string, int>::const_iterator m_iter; 
  int m_ndiff_triggers;       // Total number of HLT paths


  int m_SSlept;
};

#endif 
