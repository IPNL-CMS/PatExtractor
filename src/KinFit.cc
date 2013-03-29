/* **********************************************************************
 *                                                                      *
 * KinFit : kinematical fit for semilept tt bar events             *
 *               based on                                               *
 *               - an external  full jet parametrisation                *
 *               - chi2 minisation                                      *
 *               - computed neutrino Pz (for initialisation)            *
 *               - 2 mtop + mw constraint                               *
 * Djamel BOUMEDIENE   Djamel.Boumedien@cern.ch                         *
 *                                                                      *
 * courtesy for help and code structure :                               *
 *                           Diane Cinca diane.cinca@cern.ch            *
 * TODO : parametrization of Neutrino + Muon + electron from .dat file  *
 * internal tracking :                                                  *
 * $Id: KinFit.cc,v 1.1 2011/11/10 15:09:01 sviret Exp $*
 ************************************************************************/

#include "../interface/KinFit.h"
#include "../interface/AnalysisSettings.h"

//#define HIGH_MASS_METHOD

using namespace std;

static void * vdummy = 0;
static void MinuitFunc(int &npar, double *ddummy, double &cf, double *par, int iflag)
{
  KinFit * mySelf = (KinFit*) vdummy;
  mySelf->FuncChi2( npar, cf, par, iflag ); 
}


KinFit::KinFit(const std::string& ParamsFile, const edm::ParameterSet& cmsswSettings)
{ 
  /// argument = parametrisation file
  /// initialisations : KFChi2 kin fit chi2 set to 10^10
  /// debug level 1 by default

  m_w = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("w_mass");
  m_top = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("top_mass");
  m_b = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("b_mass");

  SigmMW = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("w_mass_error");
  SigmMT = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("top_mass_error");
  
  // Chi2 optimal values and errors
  chi2_hadronic_top_mass = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("hadronic_top_mass");

  chi2_hadronic_top_mass = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("hadronic_top_mass");
  chi2_leptonic_top_mass_semimu = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("leptonic_top_mass_semimu");
  chi2_leptonic_top_mass_semie = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("leptonic_top_mass_semie");
  chi2_hadronic_w_mass = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("hadronic_w_mass");
  chi2_pt_ttbar_system = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("pt_ttbar_system");
  chi2_ht_frac = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("ht_frac");

  chi2_sigma_hadronic_top_mass = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("sigma_hadronic_top_mass");
  chi2_sigma_leptonic_top_mass_semimu = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("sigma_leptonic_top_mass_semimu");
  chi2_sigma_leptonic_top_mass_semie = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("sigma_leptonic_top_mass_semie");
  chi2_sigma_hadronic_w_mass = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("sigma_hadronic_w_mass");
  chi2_sigma_pt_ttbar_system = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("sigma_pt_ttbar_system");
  chi2_sigma_ht_frac = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<double>("sigma_ht_frac");

  m_usePtSystInChi2 = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<bool>("use_pt_syst");
  m_useHtFracInChi2 = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<bool>("use_ht_frac");

  chi2_sigma_hadronic_top_mass_square = chi2_sigma_hadronic_top_mass * chi2_sigma_hadronic_top_mass;
  chi2_sigma_leptonic_top_mass_semimu_square = chi2_sigma_leptonic_top_mass_semimu * chi2_sigma_leptonic_top_mass_semimu;
  chi2_sigma_leptonic_top_mass_semie_square = chi2_sigma_leptonic_top_mass_semie * chi2_sigma_leptonic_top_mass_semie;
  chi2_sigma_hadronic_w_mass_square = chi2_sigma_hadronic_w_mass * chi2_sigma_hadronic_w_mass;
  chi2_sigma_pt_ttbar_system_square = chi2_sigma_pt_ttbar_system * chi2_sigma_pt_ttbar_system;
  chi2_sigma_ht_frac_square = chi2_sigma_ht_frac * chi2_sigma_ht_frac;

  if(vdummy != 0) vdummy = 0;
  KFChi2      = 10.E10;
  DEBUG_Level = -1;


  MyMinuit = new TMinuit(ParamNber);
  NStudiedFlavor=2;
  StudiedFlavor[0]=1;
  StudiedFlavor[1]=3;

  /// reads the errors from the .dat file => create one object KinFit for each type of fit
  /// This is the remaining hardcoded part which has to disappear at some point

  SigmPxNu=14.35;
  SigmPyNu=15.79;
  SigmPzNu=21.65;  	  
  SigmPTS=33.;	

  A_mu = 0.00012;
  B_mu = 0.0092;
  C_mu = 2.1;

  A_ele = 0.;
  B_ele = 0.041;
  C_ele = 1.6;

  ReadErrors(ParamsFile);
}


KinFit::~KinFit()
{
  vdummy = 0;
  delete MyMinuit;
}


int KinFit::ReadErrors(TString ParamsFile)
{
  int res=0;

  std::ifstream ofile(ParamsFile);

  if (!ofile) return 1; // file not found

  std::string var="%%";


  int security=0;

  while(var!="%%Parametrisation_starts_here%%")
  {
    ofile>>var;
    security++;
    if (security>500) return 5;
  }

  int ivar=0, lept=0;
  int ieta=0, iflavor=0; 
  double A=0, B=0, C=0, D=0;
  int ndata=0;

  while (ndata<48)
  {
    ofile>>lept>>ivar>>ieta>>iflavor>>A>>B>>C>>D;
    ndata++;

    bool ivfound=false;
    for(int iv=0; iv<NStudiedFlavor; ++iv)
    {
      if(StudiedFlavor[iv]==iflavor)
      {
        iflavor=iv;
        ivfound=true;
      }
    }

    if(!ivfound)
    {
      res=2;
      continue;
    }

    if(ieta>MaxEtaBins)
    {
      res=3;
      continue;
    }

    Errors[lept][ivar][iflavor][ieta][0]=A;
    Errors[lept][ivar][iflavor][ieta][1]=B;
    Errors[lept][ivar][iflavor][ieta][2]=C;
    Errors[lept][ivar][iflavor][ieta][3]=D;
    if ((A==0)&&(B==0)&&(C==0)&&(D==0)) res=4;

  }

  ofile.close();
  return res;
}

// ------------------------------------------------------------------------------------------
bool KinFit::ReadObjects(const TLorentzVector& Jet1, const TLorentzVector& Jet2, const TLorentzVector& BJetH, const TLorentzVector& Lepton, const TLorentzVector& Neutrino, const TLorentzVector& BJetL, bool doSemiMu, bool* eventCorrected/* = NULL*/)
{

  /*(Lepton->M()>0.05) // This is a muon    
    ? m_isMuon = true
    : m_isMuon = false;*/
  m_isMuon = doSemiMu;


  /// Measured objects are stored in  MeasuredLepton, MeasuredNeutrino, MeasuredBJetH,
  /// MeasuredBJetL, MeasuredJet1, MeasuredJet2
  MeasuredLepton   = Lepton;
  MeasuredBJetH    = BJetH;
  MeasuredBJetL    = BJetL;
  MeasuredJet1     = Jet1;
  MeasuredJet2     = Jet2;
  MeasuredNeutrino = Neutrino;


  /// We compute Pz of the neutrino

  if (! PzNeutrino(Lepton, MeasuredNeutrino, BJetL, eventCorrected))
    return false; // Don't go further if things are bad

  FittedLepton.reset(new TLorentzVector(Lepton));
  FittedBJetH.reset(new TLorentzVector(BJetH));
  FittedBJetL.reset(new TLorentzVector(BJetL));
  FittedJet1.reset(new TLorentzVector(Jet1));
  FittedJet2.reset(new TLorentzVector(Jet2));
  FittedNeutrino.reset(new TLorentzVector(MeasuredNeutrino));

  /// Sigma termes

  /// -- Utiliser cette parametrisation

  /// Resol en energie
  SigmEJet1    = DKFJetResol(Jet1, JetFlavor::W, Parameter::Energy);
  SigmEJet2    = DKFJetResol(Jet2, JetFlavor::W, Parameter::Energy);
  SigmEBJetH   = DKFJetResol(BJetH, JetFlavor::B, Parameter::Energy);
  SigmEBJetL   = DKFJetResol(BJetL, JetFlavor::B, Parameter::Energy);

  /// Resol en phi
  SigmPhiJet1  = DKFJetResol(Jet1,  JetFlavor::W, Parameter::Phi);
  SigmPhiJet2  = DKFJetResol(Jet2,  JetFlavor::W, Parameter::Phi);
  SigmPhiBJetH = DKFJetResol(BJetH, JetFlavor::B, Parameter::Phi);
  SigmPhiBJetL = DKFJetResol(BJetL, JetFlavor::B, Parameter::Phi);

  /// Resol en Eta
  SigmEtaJet1  = DKFJetResol(Jet1,  JetFlavor::W, Parameter::Eta);
  SigmEtaJet2  = DKFJetResol(Jet2,  JetFlavor::W, Parameter::Eta);
  SigmEtaBJetH = DKFJetResol(BJetH, JetFlavor::B, Parameter::Eta);
  SigmEtaBJetL = DKFJetResol(BJetL, JetFlavor::B, Parameter::Eta);

  (m_isMuon) // This is a muon
    ? SigmEMu = C_mu+B_mu*Lepton.E()+A_mu*Lepton.E()*Lepton.E()
    : SigmEMu = C_ele+B_ele*Lepton.E()+A_ele*Lepton.E()*Lepton.E();

  return true;

}

bool KinFit::Fit()
{

  if(vdummy != 0) vdummy = 0;
  vdummy = (void*) this;

  MyMinuit->mnrset(1); // Resets function value and errors to UNDEFINED

  MyMinuit->SetPrintLevel(DEBUG_Level); // debug level

  //MyMinuit->SetPrintLevel(1); // debug level
  MyMinuit->SetFCN(&MinuitFunc); // fonction to be minimised

  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1.;

  MyMinuit->mnexcm("SET ERR", arglist ,1,ierflg); 

  /// Initialize Minuit Parameters

  TString Name[ParamNber] = {
    "FitEtaJ1", "FitEtaJ2", "FitEtaBH", "FitEtaBL",
    "FitPhiJ1", "FitPhiJ2", "FitPhiBH", "FitPhiBL",
    "FitEJ1",   "FitEJ2",   "FitEBH",   "FitEBL",
    "FitEMU",
    "FitPxNu",  "FitPyNu",  "FitPzNu"};

  double vstart[ParamNber] = {
    MeasuredJet1.Eta(),    MeasuredJet2.Eta(),    MeasuredBJetH.Eta(), MeasuredBJetL.Eta(),
    MeasuredJet1.Phi(),    MeasuredJet2.Phi(),    MeasuredBJetH.Phi(), MeasuredBJetL.Phi(),
    MeasuredJet1.E(),      MeasuredJet2.E(),      MeasuredBJetH.E(),   MeasuredBJetL.E(),
    MeasuredLepton.E(),
    MeasuredNeutrino.Px(), MeasuredNeutrino.Py(), MeasuredNeutrino.Pz()
  };

  double step[ParamNber];
  double LowerBound[ParamNber];
  double UpperBound[ParamNber];

  for(int i=0 ; i < ParamNber ; ++i )
  {
    (i < 8)
      ? step[i]=1E-6
      : step[i]=1E-4;

    if (i < 4) // Eta bound
    {
      LowerBound[i]= vstart[i] - 2.;
      UpperBound[i]= vstart[i] + 2.;
    }    

    if (i >= 4 && i < 8) // Phi bound
    {
      LowerBound[i] = vstart[i] - 1.;
      UpperBound[i] = vstart[i] + 1.;
    }

    if (i >= 8 && i < 13) // NRJ bound
    {
      LowerBound[i] = 0.5 * vstart[i];
      UpperBound[i] = 1.5 * vstart[i];
    }

    if (i >= 13 && i < 16) // Neutrino p bound
    {
      LowerBound[i] = vstart[i] - 0.5 * fabs(vstart[i]);
      UpperBound[i] = vstart[i] + 0.5 * fabs(vstart[i]);
    }

    MyMinuit->mnparm(i, Name[i], vstart[i], step[i], LowerBound[i], UpperBound[i], ierflg);
    if(ierflg)
    {
      std::cout << "PROBLEM WITH PARAMETER " << i << "\t => ierflg=" << ierflg << std::endl;
      std::cout << "ParamNum= " << i << "\tName = " << Name[i] 
        << "\tvstart=" << vstart[i] << "\tstep="<< step[i] 
        << std::endl;
      break;
    }
  }

  arglist[0] = 1000.;
  arglist[1] = 1.;

  /// minimisation
  MyMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  /// Print results
  vdummy = 0;
  double amin,edm,errdef;
  int nvpar,nparx,icstat;

  MyMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  /// Get parameters
  for(int i = 0; i < ParamNber; ++i)
  {
    MyMinuit->GetParameter(i, FittedParam[i], ErrFitParam[i]);
  }

  FuncChi2(ParamNber, KFChi2, FittedParam, 0); 

  // status = 1 if minuit output is "CONVERGED";
  return (MyMinuit->fCstatu == "CONVERGED ");
}


// Change this, and use LorentzVector from ROOT instead
// See: http://project-mathlibs.web.cern.ch/project-mathlibs/sw/html/Vector.html


/// Function chi2 dont le pointeur est fourni a minuit pour minimisation
void KinFit::FuncChi2(const int &npar, double &f, double *par, int iflag)
{
  FittedJet1.reset(new TLorentzVector(par[8]*TMath::Cos(par[4])/TMath::CosH(par[0]), par[8]*TMath::Sin(par[4])/TMath::CosH(par[0]), par[8]*TMath::TanH(par[0]), par[8]));

  FittedJet2.reset(new TLorentzVector(par[9]*TMath::Cos(par[5])/TMath::CosH(par[1]), par[9]*TMath::Sin(par[5])/TMath::CosH(par[1]), par[9]*TMath::TanH(par[1]),par[9]));

  FittedBJetH.reset(new TLorentzVector(sqrt((par[10]*par[10]-m_b*m_b))*TMath::Cos(par[6])/TMath::CosH(par[2]), sqrt((par[10]*par[10]-m_b*m_b))*TMath::Sin(par[6])/TMath::CosH(par[2]), sqrt((par[10]*par[10]-m_b*m_b))*TMath::TanH(par[2]),par[10]));

  FittedBJetL.reset(new TLorentzVector(sqrt((par[11]*par[11]-m_b*m_b))*TMath::Cos(par[7])/TMath::CosH(par[3]), sqrt((par[11]*par[11]-m_b*m_b))*TMath::Sin(par[7])/TMath::CosH(par[3]), sqrt((par[11]*par[11]-m_b*m_b))*TMath::TanH(par[3]),par[11]));

  FittedLepton.reset(new TLorentzVector(par[12]*TMath::Cos(MeasuredLepton.Phi())/TMath::CosH(MeasuredLepton.Eta()), par[12]*TMath::Sin(MeasuredLepton.Phi())/TMath::CosH(MeasuredLepton.Eta()), par[12]*TMath::TanH(MeasuredLepton.Eta()),par[12]));

  FittedNeutrino.reset(new TLorentzVector(par[13], par[14], par[15], sqrt(par[13] * par[13] + par[14] * par[14] + par[15] * par[15])));
  
  f = Chi2(); 
}


double KinFit::Chi2()
{
  double ResJet1Eta    = (*FittedJet1).Eta()-(MeasuredJet1).Eta();
  double ResJet2Eta    = (*FittedJet2).Eta()-(MeasuredJet2).Eta();
  double ResBJetHEta   = (*FittedBJetH).Eta()-(MeasuredBJetH).Eta();
  double ResBJetLEta   = (*FittedBJetL).Eta()-(MeasuredBJetL).Eta();
  double ResJet1Phi    = (*FittedJet1).Phi()-(MeasuredJet1).Phi();
  double ResJet2Phi    = (*FittedJet2).Phi()-(MeasuredJet2).Phi();
  double ResBJetHPhi   = (*FittedBJetH).Phi()-(MeasuredBJetH).Phi();
  double ResBJetLPhi   = (*FittedBJetL).Phi()-(MeasuredBJetL).Phi();
  double ResJet1E      = (*FittedJet1).E()-(MeasuredJet1).E();
  double ResJet2E      = (*FittedJet2).E()-(MeasuredJet2).E();
  double ResBJetHE     = (*FittedBJetH).E()-(MeasuredBJetH).E();
  double ResBJetLE     = (*FittedBJetL).E()-(MeasuredBJetL).E();   
  double ResLeptonE    = (*FittedLepton).E()-(MeasuredLepton).E();
  double ResNeutrinoPx = (*FittedNeutrino).Px()-(MeasuredNeutrino).Px();
  double ResNeutrinoPy = (*FittedNeutrino).Py()-(MeasuredNeutrino).Py();
  double ResNeutrinoPz = (*FittedNeutrino).Pz()-(MeasuredNeutrino).Pz();


  /// Compute Masses
  MjjViaFit   = (*FittedJet1+*FittedJet2).M();
  MlvViaFit   = (*FittedLepton+*FittedNeutrino).M();
  MjjbViaFit  = (*FittedJet1+*FittedJet2+*FittedBJetH).M();
  MlvbViaFit  = (*FittedLepton+*FittedNeutrino+*FittedBJetL).M();
  FitPtSystem = (*FittedJet1+*FittedJet2+*FittedBJetH+*FittedLepton+*FittedNeutrino+*FittedBJetL).Pt();

  /// Compute Chi2
  double finalchi=0.;

  finalchi+=TMath::Power(ResJet1Eta/SigmEtaJet1,2.);
  finalchi+=TMath::Power(ResJet2Eta/SigmEtaJet2,2.);
  finalchi+=TMath::Power(ResBJetHEta/SigmEtaBJetH,2.);
  finalchi+=TMath::Power(ResBJetLEta/SigmEtaBJetL,2.);
  finalchi+=TMath::Power(ResJet1Phi/SigmPhiJet1,2.);
  finalchi+=TMath::Power(ResJet2Phi/SigmPhiJet2,2.);
  finalchi+=TMath::Power(ResBJetHPhi/SigmPhiBJetH,2.);
  finalchi+=TMath::Power(ResBJetLPhi/SigmPhiBJetL,2.);
  finalchi+=TMath::Power(ResJet1E/SigmEJet1,2.);
  finalchi+=TMath::Power(ResJet2E/SigmEJet2,2.);
  finalchi+=TMath::Power(ResBJetHE/SigmEBJetH,2.);
  finalchi+=TMath::Power(ResBJetLE/SigmEBJetL,2.);
  finalchi+=TMath::Power(ResLeptonE/SigmEMu,2.);
  finalchi+=TMath::Power(ResNeutrinoPx/SigmPxNu,2.);
  finalchi+=TMath::Power(ResNeutrinoPy/SigmPyNu,2.);
  finalchi+=TMath::Power(ResNeutrinoPz/SigmPzNu,2.);
  finalchi+=TMath::Power((MjjViaFit-m_w)/SigmMW,2.);
  finalchi+=TMath::Power((MlvViaFit-m_w)/SigmMW,2.);
  finalchi+=TMath::Power((MjjbViaFit-m_top)/SigmMT,2.);
  finalchi+=TMath::Power((MlvbViaFit-m_top)/SigmMT,2.); 

  return finalchi;
}

//viola: i duplicate float Chi2::GlobalSimpleChi2 to have a second implementation of the method 
//that modifies a tlorentzvector for met
double KinFit::GlobalSimpleChi2(double totPt)
{


  /*
  std::cout << "MET: " << MeasuredNeutrino.Pt() << "; " << MeasuredNeutrino.Px() << "; " << MeasuredNeutrino.Py() << "; " << MeasuredNeutrino.Pz() << "; " << MeasuredNeutrino.E() << std::endl;
  std::cout << "Lepton: " << (*MeasuredLepton).Pt() << std::endl;

  std::cout << "Jet 1: " << (*MeasuredJet1).Pt() << std::endl;
  std::cout << "Jet 2: " << (*MeasuredJet2).Pt() << std::endl;

  std::cout << "Jet B-H: " << (*MeasuredBJetH).Pt() << std::endl;
  std::cout << "Jet B-L: " << (*MeasuredBJetL).Pt() << std::endl;

  std::cout << "Total pt: " << totPt << std::endl;
  */

  float MW          = sqrt(max(0., (MeasuredJet1 + MeasuredJet2).M2()));
  //std::cout << "MW: " << MW << std::endl;
  float MtopH       = sqrt(max(0., (MeasuredJet1 + MeasuredJet2 + MeasuredBJetH).M2()));
  //std::cout << "MTopH: " << MtopH << std::endl;
  float MtopL       = sqrt(max(0., (MeasuredNeutrino + MeasuredLepton + MeasuredBJetL).M2()));
  //std::cout << "MTopL: " << MtopL << std::endl;
  //float SolPtSystem = (*MeasuredJet1+*MeasuredJet2+*MeasuredBJetL+*MeasuredBJetH).Pt()/totPt;
  float SolPtSystem = (MeasuredJet1.Pt() + MeasuredJet2.Pt() + MeasuredBJetL.Pt() + MeasuredBJetH.Pt()) / totPt;
  //std::cout << "Pt syst: " << SolPtSystem << std::endl;
  float TTbarSystemPt = ((MeasuredJet1 + MeasuredJet2 + MeasuredBJetL + MeasuredBJetH + MeasuredNeutrino + MeasuredLepton).Pt());

  float chi2 = ((MtopH - chi2_hadronic_top_mass) * (MtopH - chi2_hadronic_top_mass) / (chi2_sigma_hadronic_top_mass_square)) + 
    ((MW - chi2_hadronic_w_mass) * (MW - chi2_hadronic_w_mass) / (chi2_sigma_hadronic_w_mass_square));

  if (m_usePtSystInChi2) {
    chi2 += ((SolPtSystem - chi2_ht_frac) * (SolPtSystem - chi2_ht_frac) / (chi2_sigma_ht_frac_square));
  }

  if (m_useHtFracInChi2) {
    chi2 += ((TTbarSystemPt - chi2_pt_ttbar_system) * (TTbarSystemPt - chi2_pt_ttbar_system) / chi2_sigma_pt_ttbar_system_square);
  }

  /// lepton dependant part
  (m_isMuon)
    ? chi2 += (MtopL - chi2_leptonic_top_mass_semimu) * (MtopL - chi2_leptonic_top_mass_semimu) / (chi2_sigma_leptonic_top_mass_semimu_square)
    : chi2 += (MtopL - chi2_leptonic_top_mass_semie) * (MtopL - chi2_leptonic_top_mass_semie) / (chi2_sigma_leptonic_top_mass_semie_square);

  return  chi2;
}

double KinFit::GetKFChi2()
{
  return KFChi2;
}

const int KinFit::GetParamNber()
{
  return ParamNber;
}

double KinFit::GetFitVar(const int Num)
{
  return FittedParam[Num];
}

double KinFit::GetErrorFitVar(const int Num)
{
  return ErrFitParam[Num];
}



// -----------------------------------------------------------------------------------------------//
// DKFJetResol provides an uncertainty for each jet                                               //
// input : TLorentzVector                                                                         //
//         JetFlavor : 1 or 2 (W or b)                                                            //
//         IPar : parameter for which the uncertainty is computed                                 //
//              = 1, 2 or 3 (Energy, Eta or Phi)                                                  //
// output : Absolute Uncertainty on the parameter IPar                                            //
// Remark : This function requires Errors[][][][] to be filled                               //
//          by the function ReadErrors() called once when the                                //
//          object KinFit is created                                                         //
// -----------------------------------------------------------------------------------------------//

// arg1 : TLorentzVector, arg2: jet type (1=W, 2=b), arg3: type de resolution 1=Ene, 2=Eta, 3=Phi	   

double KinFit::DKFJetResol(const TLorentzVector& Jet, JetFlavor flavor, Parameter IPar)
{
  double res=1E10;

  int ind_lept = 0;

  if (!m_isMuon) ind_lept = 1; // Error index for electrons

/*
  if(IPar>3)
  {
    std::cout<<" <!> ERROR in DKFJetResol. Unknown parameter "<<IPar<<". Parameter number is limitted to 1 -> 3"<<std::endl;
    return res;
  }
*/

  if (fabs(Errors[ind_lept][0][0][0][0])<1E-10 &&
      fabs(Errors[ind_lept][0][0][0][1])<1E-10 &&
      fabs(Errors[ind_lept][0][0][0][2])<1E-10 &&
      fabs(Errors[ind_lept][0][0][0][3])<1E-10)
  {
    /// A=B=C=D=0 for E parametrisation (Wj) in all Eta bins
    std::cout<<" <!> ERROR in KinFit. Errors are not loaded ! Error set to 1E10 ..."<<std::endl;
    std::cout<<" <!> Check if the text file containing the errors is there"<<std::endl;
    return res;
  }

  /// parametres pour JetW
  double A=0, B=0, C=0, D=0;

  double Eta = Jet.Eta();
  double Ene = Jet.E();

  unsigned int ieta=0;


  /// barrel
  if(fabs(Eta)<1.4) ieta = 0; 

  /// crak
  else if((fabs(Eta)>=1.4)&&(fabs(Eta)<1.7)) ieta = 1; 

  /// end cap
  else if((fabs(Eta)>=1.7)&&(fabs(Eta)<2.6)) ieta = 2; 

  /// tres bas angle
  else if(fabs(Eta)>=2.6) ieta = 3;

  /// nouvelles (par fichier .dat)
  A=Errors[ind_lept][static_cast<int>(IPar)][static_cast<int>(flavor)][ieta][0];
  B=Errors[ind_lept][static_cast<int>(IPar)][static_cast<int>(flavor)][ieta][1];
  C=Errors[ind_lept][static_cast<int>(IPar)][static_cast<int>(flavor)][ieta][2];
  D=Errors[ind_lept][static_cast<int>(IPar)][static_cast<int>(flavor)][ieta][3];

  (Ene>0)
    ? res= A+B*sqrt(Ene)+C/Ene+D/(Ene*Ene)
    : res= A;

  return res;
}


TLorentzVector* KinFit::GetFittedFirstLightJet()
{
  return FittedJet1.get();
}

TLorentzVector* KinFit::GetFittedSecondLightJet()
{
  return FittedJet2.get();
}

TLorentzVector* KinFit::GetFittedLeptonicBJet()
{
  return FittedBJetL.get();
}

TLorentzVector* KinFit::GetFittedHadronicBJet()
{
  return FittedBJetH.get();
}

// Backward compatibility
TLorentzVector* KinFit::GetFittedJet(int index)
{
  switch ( index ) {
    case 0: return GetFittedFirstLightJet();
    case 1: return GetFittedSecondLightJet();
    case 2: return GetFittedHadronicBJet();
    case 3: return GetFittedLeptonicBJet();
  }

  return nullptr;
}


TLorentzVector* KinFit::GetFittedLepton()
{
  return FittedLepton.get();
}



TLorentzVector* KinFit::GetFittedNeutrino()
{
  return FittedNeutrino.get();
}


///////////////////////////////
// Quelques fonctions utiles //
///////////////////////////////


double KinFit::PzNeutrino(const TLorentzVector& lept, TLorentzVector& neut, const TLorentzVector& bJet, bool* eventCorrected/* = NULL*/)
{
  if (eventCorrected)
    *eventCorrected = false;

  if(!lept.E()) return 0;

  double x = (m_w * m_w - lept.M() * lept.M() + 2. * (neut.Px() * lept.Px() + neut.Py() * lept.Py())) / (2 * lept.E());
  double a = 1 - (lept.Pz() * lept.Pz()) / (lept.E() * lept.E());
  double b = -2. * (lept.Pz() / lept.E()) * x;
  double c = neut.Pt() * neut.Pt() - x * x;

  if (!a && !b) return 0;

  if(!a)
  {     
    neut.SetPz(-1 * c / b);
    neut.SetE(sqrt(neut.Px() * neut.Px() + neut.Py() * neut.Py() + neut.Pz() * neut.Pz()));
    return 1;
  }


  double delta = b * b - 4 * a *c;

#ifdef HIGH_MASS_METHOD
  if (delta < 0) {
    if (eventCorrected)
      *eventCorrected = true;

    delta = 0;
  }
#else
  if (delta < 0)     // No solution, try to correct MET
  {
    double rat = neut.Py() / neut.Px();

    double u = 4./(lept.E()*lept.E())*((lept.Px()+rat*lept.Py())*(lept.Px()+rat*lept.Py())/(1+rat*rat)
        -(lept.E()*lept.E())+(lept.Pz()*lept.Pz()));

    double v = 4./(lept.E()*lept.E())*(m_w*m_w-lept.M()*lept.M())
      *(lept.Px()+rat*lept.Py())/sqrt(1+rat*rat);

    double w = (m_w*m_w-lept.M()*lept.M())*(m_w*m_w-lept.M()*lept.M())/(lept.E()*lept.E());

    double deltan = v * v - 4 * u * w;

    if (deltan < 0)
      return 0; // Hopeless, MET can't be corrected

    double pt      = 0.;
    double corfact = 0.;

    if(u==0)
    {
      pt = -w/v;
      if (pt <= 0) return 0; // There is no way out...

      corfact = pt / neut.Pt();
    }
    else // Deltan>=0 and u!=0
    {
      double pt2 = (v-(sqrt(deltan)))/(2*u);
      double pt1 = (v+(sqrt(deltan)))/(2*u);

      // Pas de correction car negative
      if(pt1<=0 && pt2<=0) return 0;

      if(pt1>0 && pt2<0) pt=pt1;
      if(pt2>0 && pt1<0) pt=pt2;
      if(pt1>0 && pt2>0)
      {
        (fabs(pt1-neut.Pt())<=fabs(pt2-neut.Pt()))
          ? pt=pt1
          : pt=pt2;     
      }

      corfact = pt/neut.Pt();
    }

    // Now we have the correction factor

    neut.SetPx(corfact*neut.Px());
    neut.SetPy(corfact*neut.Py());

    // Recompute the new parameters

    x = (m_w*m_w-lept.M()*lept.M()+2.*(neut.Px()*lept.Px()+neut.Py()*lept.Py()))/(2*lept.E());
    a = 1-(lept.Pz()*lept.Pz())/(lept.E()*lept.E());
    b = -2.*(lept.Pz()/lept.E())*x;
    c = neut.Px()*neut.Px()+neut.Py()*neut.Py()-x*x;

//         std::cout << "We have rescaled the MET " << lept->E() << " / " << corfact <<  " , now delta should be null:" << std::endl;
//         std::cout << "Previous delta: " << delta<< std::endl;


    delta= b*b-4*a*c;

    if (fabs(delta)<0.000001) delta=0.;

//    std::cout << "New delta     : " << delta << std::endl;

    if (delta != 0)
      return 0; // This should not happen, but who knows...

    if (eventCorrected)
      *eventCorrected = true;
  }
#endif


  // We can go back to the normal path: 

  TLorentzVector TopCand1 = lept + bJet;
  TLorentzVector TopCand2 = lept + bJet;

  neut.SetPz((-b-(sqrt(delta)))/(2*a));
  neut.SetE(sqrt(neut.Px()*neut.Px()+neut.Py()*neut.Py()+neut.Pz()*neut.Pz()));
  TopCand1 += neut;

  neut.SetPz((-b+(sqrt(delta)))/(2*a));
  neut.SetE(sqrt(neut.Px()*neut.Px()+neut.Py()*neut.Py()+neut.Pz()*neut.Pz()));
  TopCand2 += neut;

  double mtt_1 = sqrt(std::max(0.,TopCand1.M2()));
  double mtt_2 = sqrt(std::max(0.,TopCand2.M2()));

  if(fabs(mtt_1-m_top) <= fabs(mtt_2-m_top)) // Otherwise it's OK
  {
    neut.SetPz((-b-(sqrt(delta)))/(2*a));
    neut.SetE(sqrt(neut.Px()*neut.Px()+neut.Py()*neut.Py()+neut.Pz()*neut.Pz()));
  }

  return 1;
} 

