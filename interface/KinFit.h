#ifndef KINFIT_H
#define KINFIT_H

/// ROOT includes
#include "TLorentzVector.h"
#include "TMinuit.h"
#include <TROOT.h>
#include <TMath.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>

#include <memory> // for std::shared_ptr

/// modif BDE, reduction nombre de parametre car mtop fixee
static const int ParamNber = 16;
static const int MaxEtaBins=4;



class KinFit
{ 
public:
  

  KinFit(TString ParamsFile
	 , double wmass, double topmass, double bmass
	 , double wmass_e, double topmass_e); 
  virtual ~KinFit();
  
  double Chi2();
  double GlobalSimpleChi2(double TotPt);
  
  int ReadObjects(TLorentzVector* Jet1, TLorentzVector* Jet2, TLorentzVector* BJetH, TLorentzVector* Lepton, TLorentzVector* Neutrino, TLorentzVector* BJetL, bool doSemiMu);
  
  bool Fit();
  
  double GetKFChi2();
  double GetFitVar(const int Num);
  double GetErrorFitVar(const int Num);
  const int GetParamNber();
  inline void SetDebugMode(int DBL){ DEBUG_Level = DBL; return; }
  int ReadErrors(TString ParamsFile);
  double DKFJetResol(TLorentzVector *Jet,int JetFlavor, int IPar);
  
  TLorentzVector* GetFittedFirstLightJet();
  TLorentzVector* GetFittedSecondLightJet();
  TLorentzVector* GetFittedLeptonicBJet();
  TLorentzVector* GetFittedHadronicBJet();
  TLorentzVector* GetFittedJet(int index);
  TLorentzVector* GetFittedLepton();
  TLorentzVector* GetFittedNeutrino();
  
  double PzNeutrino(TLorentzVector *lept, TLorentzVector *neut, TLorentzVector *bJet);
  
  void FuncChi2(const int &npar, double &f, double *par, int iflag);
  
 private:    
  TMinuit *MyMinuit;
  
  bool m_isMuon;
  
  //const measured termes : input in Fit
  double KFChi2;
  
  double MeanMeasMTop;
  
  /// Computed termes
  double MjjViaFit;
  double MlvViaFit;
  double MjjbViaFit;
  double MlvbViaFit;
  
  /// Other Needed 
  double FitPtSystem;
  
  double FittedParam[ParamNber];
  double ErrFitParam[ParamNber];
  
  /// Sigma termes
  double SigmEtaJet1;
  double SigmEtaJet2;
  double SigmEtaBJetH;
  double SigmEtaBJetL;
  double SigmPhiJet1;
  double SigmPhiJet2;
  double SigmPhiBJetH;
  double SigmPhiBJetL;
  double SigmEJet1;
  double SigmEJet2;
  double SigmEBJetH;
  double SigmEBJetL;
  double SigmEMu;
  
  double SigmPxNu;
  double SigmPyNu;
  double SigmPzNu;
  double SigmMW;
  double SigmMT;
  double SigmPTS;    
  double A_mu;
  double B_mu;
  double C_mu;
  double A_ele;
  double B_ele;
  double C_ele;
  
  /// constant masses
  double m_w;  
  double m_b;  
  double m_top;
  
  
  int DEBUG_Level;
  
  double Errors[2][3][3][4][4];
  /// This vector contains a set of parameters used to parametrise the errors.
  /// It is filled by ReadErrors() which should be called once in the constructor
  
  /// Studied flavors
  /// only two flavors are considered : Wjets=1, bjets=3
  int NStudiedFlavor;
  int StudiedFlavor[2];
  //StudiedFlavor[0]=1;
  //StudiedFlavor[1]=3;
  
  //double MeasParams[30];
  TLorentzVector MeasuredNeutrino;

  TLorentzVector* MeasuredLepton;
  TLorentzVector* MeasuredBJetH;
  TLorentzVector* MeasuredBJetL;
  TLorentzVector* MeasuredJet1;
  TLorentzVector* MeasuredJet2;
  
  std::shared_ptr<TLorentzVector> FittedLepton;
  std::shared_ptr<TLorentzVector> FittedNeutrino;
  std::shared_ptr<TLorentzVector> FittedBJetH;
  std::shared_ptr<TLorentzVector> FittedBJetL;
  std::shared_ptr<TLorentzVector> FittedJet1;
  std::shared_ptr<TLorentzVector> FittedJet2;

};


#endif

