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

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <memory> // for std::shared_ptr
#include <cstdint>

/// modif BDE, reduction nombre de parametre car mtop fixee
static const int ParamNber = 16;
static const int MaxEtaBins=4;

enum class JetFlavor: std::uint8_t {
  W = 0, B = 1
};

enum class Parameter: std::uint8_t {
  Energy = 0, Eta = 1, Phi = 2
};

class AnalysisSettings;


class KinFit
{ 
public:
  

  KinFit(const std::string& ParamsFile, const edm::ParameterSet& settings);
  virtual ~KinFit();
  
  double Chi2();
  double GlobalSimpleChi2(double TotPt);
  
  bool ReadObjects(const TLorentzVector& Jet1, const TLorentzVector& Jet2, const TLorentzVector& BJetH, const TLorentzVector& Lepton, const TLorentzVector& Neutrino, const TLorentzVector& BJetL, bool doSemiMu, bool* eventCorrected = NULL);
  
  bool Fit();
  
  double GetKFChi2();
  double GetFitVar(const int Num);
  double GetErrorFitVar(const int Num);
  const int GetParamNber();
  inline void SetDebugMode(int DBL){ DEBUG_Level = DBL; return; }
  int ReadErrors(TString ParamsFile);

  double DKFJetResol(const TLorentzVector& Jet, JetFlavor flavor, Parameter IPar);
  
  TLorentzVector* GetFittedFirstLightJet();
  TLorentzVector* GetFittedSecondLightJet();
  TLorentzVector* GetFittedLeptonicBJet();
  TLorentzVector* GetFittedHadronicBJet();
  TLorentzVector* GetFittedJet(int index);
  TLorentzVector* GetFittedLepton();
  TLorentzVector* GetFittedNeutrino();

  const TLorentzVector& GetMeasuredFirstLightJet() {
    return MeasuredJet1;
  }

  const TLorentzVector& GetMeasuredSecondLightJet() {
    return MeasuredJet2;
  }

  const TLorentzVector& GetMeasuredLeptonicBJet() {
    return MeasuredBJetL;
  }

  const TLorentzVector& GetMeasuredHadronicBJet() {
    return MeasuredBJetH;
  }

  const TLorentzVector& GetMeasuredLepton() {
    return MeasuredLepton;
  }

  const TLorentzVector& GetMeasuredNeutrino() {
    return MeasuredNeutrino;
  }
  
  double PzNeutrino(const TLorentzVector& lept, TLorentzVector& neut, const TLorentzVector& bJet, bool* eventCorrected = NULL);
  
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
  
  // Chi2
  double chi2_hadronic_top_mass;
  double chi2_leptonic_top_mass_semimu;
  double chi2_leptonic_top_mass_semie;
  double chi2_hadronic_w_mass;
  double chi2_pt_ttbar_system;
  double chi2_ht_frac;

  double chi2_sigma_hadronic_top_mass;
  double chi2_sigma_leptonic_top_mass_semimu;
  double chi2_sigma_leptonic_top_mass_semie;
  double chi2_sigma_hadronic_w_mass;
  double chi2_sigma_pt_ttbar_system;
  double chi2_sigma_ht_frac;

  double chi2_sigma_hadronic_top_mass_square;
  double chi2_sigma_leptonic_top_mass_semimu_square;
  double chi2_sigma_leptonic_top_mass_semie_square;
  double chi2_sigma_hadronic_w_mass_square;
  double chi2_sigma_pt_ttbar_system_square;
  double chi2_sigma_ht_frac_square;
  
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
  TLorentzVector MeasuredLepton;
  TLorentzVector MeasuredBJetH;
  TLorentzVector MeasuredBJetL;
  TLorentzVector MeasuredJet1;
  TLorentzVector MeasuredJet2;
  
  std::shared_ptr<TLorentzVector> FittedLepton;
  std::shared_ptr<TLorentzVector> FittedNeutrino;
  std::shared_ptr<TLorentzVector> FittedBJetH;
  std::shared_ptr<TLorentzVector> FittedBJetL;
  std::shared_ptr<TLorentzVector> FittedJet1;
  std::shared_ptr<TLorentzVector> FittedJet2;

  // Chi2 testing
  bool m_usePtSystInChi2;
  bool m_useHtFracInChi2;
};


#endif

