#ifndef DIMUON_ANALYSIS_H
#define DIMUON_ANALYSIS_H

/*****************************

Simple example class showing how to perform an 
analysis using the PatExtractor tools

S.Viret (viret@in2p3.fr): 31/05/11

More info: http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.PHYTuto

 ******************************/

//Include std C++
#include <iostream>
#include <vector>
using namespace std;

#include "AnalysisSettings.h"
#include "MuonExtractor.h"
#include "TLorentzVector.h"

class dimuon_analysis
{
 public:
  dimuon_analysis(AnalysisSettings *settings);

  ~dimuon_analysis();
  
  //Selection

  bool isMuSel(TLorentzVector *muon);
  int  dimuon_Sel(MuonExtractor *muon, int evtnum);


  void reset();
  void fillTree();
   
 private:

  TTree* m_tree_dimuon;

  int m_evt;

  float m_mdimuon;
  float	m_totp;
  float	m_mu1pt;
  float	m_mu2pt;

  float	m_mu_Ptcut;
  float	m_mu_Mult;
};

#endif 
