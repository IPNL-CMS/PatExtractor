#include "SingleTprime_analysis.h"

#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

#include "TH2.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "Extractors/PatExtractor/interface/AnalysisSettings.h"
#include "Extractors/PatExtractor/interface/MCExtractor.h"
#include "Extractors/PatExtractor/interface/HLTExtractor.h"
#include "Extractors/PatExtractor/interface/MuonExtractor.h"
#include "Extractors/PatExtractor/interface/ElectronExtractor.h"
#include "Extractors/PatExtractor/interface/METExtractor.h"
#include "Extractors/PatExtractor/interface/VertexExtractor.h"
#include "Extractors/PatExtractor/interface/KinFit.h"
#include "Extractors/PatExtractor/interface/EventExtractor.h"
#include "Extractors/PatExtractor/interface/PatExtractor.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Constants 
int NumberOfGoodJets=5;
int NumberOfBadJets=1;
float THTcut=630;
float HiggsMass=125.0;
float HiggsMassWindow=2000.0;
float WMass=80.3;
float WMassWindow=2000.0;
float TopMass=172.5;
float TopMassWindow=2000.0;

using namespace std;

namespace patextractor {

  SingleTprime_analysis::SingleTprime_analysis(const edm::ParameterSet& iConfig): Plugin(iConfig)/*,
    jetEnergyResolutionScaleFactors_  (cmsswSettings.getParameter<std::vector<double> >("jetEnergyResolutionScaleFactors")),
    jetEnergyResolutionEtaBinning_    (cmsswSettings.getParameter<std::vector<double> >("jetEnergyResolutionEtaBinning"))*/
{
  std::cout << "Entering SingleTprime analysis" << std::endl;

  // Lorentz vectors

  ReconstructedHiggs = new TLorentzVector(0., 0., 0., 0.);
  ReconstructedW = new TLorentzVector(0., 0., 0., 0.);
  ReconstructedTop = new TLorentzVector(0., 0., 0., 0.);
  ReconstructedTprime = new TLorentzVector(0., 0., 0., 0.);
  TrueHiggs = new TLorentzVector(0., 0., 0., 0.);
  TrueW = new TLorentzVector(0., 0., 0., 0.);
  TrueTop = new TLorentzVector(0., 0., 0., 0.);
  TrueTprime = new TLorentzVector(0., 0., 0., 0.);
  FirstHiggsJet = new TLorentzVector(0., 0., 0., 0.);
  SecondHiggsJet = new TLorentzVector(0., 0., 0., 0.);
  FirstWJet = new TLorentzVector(0., 0., 0., 0.);
  SecondWJet = new TLorentzVector(0., 0., 0., 0.);
  TopJet = new TLorentzVector(0., 0., 0., 0.);
  FirstTrueHiggsJet = new TLorentzVector(0., 0., 0., 0.);
  SecondTrueHiggsJet = new TLorentzVector(0., 0., 0., 0.);
  FirstTrueWJet = new TLorentzVector(0., 0., 0., 0.);
  SecondTrueWJet = new TLorentzVector(0., 0., 0., 0.);
  TopTrueJet = new TLorentzVector(0., 0., 0., 0.);

  /// Tree definition containing your analysis results

  m_tree_stp = new TTree("stp","Single Analysis info");  
  
  /// Branches definition

  m_tree_stp->Branch("Reconstructed_Higgs", &ReconstructedHiggs);
  m_tree_stp->Branch("Reconstructed_W", &ReconstructedW);
  m_tree_stp->Branch("Reconstructed_Top", &ReconstructedTop);
  m_tree_stp->Branch("Reconstructed_Tprime", &ReconstructedTprime);  
  m_tree_stp->Branch("True_Higgs", &TrueHiggs);
  m_tree_stp->Branch("True_W", &TrueW);
  m_tree_stp->Branch("True_Top", &TrueTop);
  m_tree_stp->Branch("True_Tprime", &TrueTprime);
  m_tree_stp->Branch("First_Higgs_Jet", &FirstHiggsJet);
  m_tree_stp->Branch("Second_Higgs_Jet", &SecondHiggsJet);
  m_tree_stp->Branch("First_W_Jet", &FirstWJet);
  m_tree_stp->Branch("Second_W_Jet", &SecondWJet);
  m_tree_stp->Branch("Top_Jet", &TopJet);
  m_tree_stp->Branch("First_True_Higgs_Jet", &FirstTrueHiggsJet);
  m_tree_stp->Branch("Second_True_Higgs_Jet", &SecondTrueHiggsJet);
  m_tree_stp->Branch("First_True_W_Jet", &FirstTrueWJet);
  m_tree_stp->Branch("Second_True_W_Jet", &SecondTrueWJet);
  m_tree_stp->Branch("True_Top_Jet", &TopTrueJet);
  m_tree_stp->Branch("evt",         &m_evt    ,"evt/I");      // Simple evt number or event ID
  m_tree_stp->Branch("PU",         &m_nPU    ,"m_nPU/I");
  m_tree_stp->Branch("THT", &m_THT,"THT/F");
  m_tree_stp->Branch("jet1_pt",  &m_jet1pt   ,"jet1_pt/F");
  m_tree_stp->Branch("jet2_pt",  &m_jet2pt   ,"jet2_pt/F");
  m_tree_stp->Branch("jet3_pt",  &m_jet3pt   ,"jet3_pt/F");
  m_tree_stp->Branch("jet4_pt",  &m_jet4pt   ,"jet4_pt/F");
  m_tree_stp->Branch("jet5_pt",  &m_jet5pt   ,"jet5_pt/F");
  m_tree_stp->Branch("jet6_pt",  &m_jet6pt   ,"jet6_pt/F");
  m_tree_stp->Branch("Number_CSVLbtagged_jets",  &m_NBtaggedJets_CSVL   ,"m_NBtaggedJets_CSVL/I");
  m_tree_stp->Branch("DeltaR_of_Higgs_Jets",  &m_DRHiggsJets   ,"DRHiggsJets/F");
  m_tree_stp->Branch("DeltaR_of_W_Jets",  &m_DRWJets   ,"DRWJets/F");
  m_tree_stp->Branch("DeltaR_of_Top_Higgs",  &m_DRTopHiggs   ,"DRTopHiggs/F");
  m_tree_stp->Branch("DeltaR_of_W_Higgs",  &m_DRWHiggs   ,"DRWHiggs/F");
  m_tree_stp->Branch("Relative_THT",  &m_RelTHT   ,"RelTHT/F");
  m_tree_stp->Branch("Correct_Tprime",  &CorrectTprime   ,"CorrectTprime/I");
  m_tree_stp->Branch("Correct_Higgs",  &CorrectH   ,"CorrectH/I");
  m_tree_stp->Branch("Correct_W",  &CorrectW   ,"CorrectW/I");
  m_tree_stp->Branch("Correct_Top",  &CorrectTop   ,"CorrectTop/I");
  m_tree_stp->Branch("Correct_Higgs_Jet",  &CorrectHiggsJet   ,"CorrectHiggsJet/I");
  m_tree_stp->Branch("Correct_W_Jet",  &CorrectWJet   ,"CorrectWJet/I");
  m_tree_stp->Branch("Correct_Top_Jet",  &CorrectTopJet   ,"CorrectTopJet/I");

  // Initialize the analysis parameters using the ParameterSet iConfig
  //int an_option = iConfig.getUntrackedParameter<int>("an_option", 0);
  m_jet_Ptcut = 30;                               // Default val
  m_jet_EtaMaxcut = 5.0;
  m_jet_EtaAccepcut = 2.5;
  m_jet_OverlapAccep = 0.5;
  m_jet_MultInAcceptance = 5;
  m_jet_MultOutAcceptance = 1;
  m_JET_btag_CSVL = 0.244;
  evt_num = 0;
  m_DRMatching=0.3;
  m_DPtMatching=10.0;
}

SingleTprime_analysis::~SingleTprime_analysis(){}

  int SingleTprime_analysis::SingleTprime_Sel() //Main function for the analysis
{
  int n_jets = m_jetMet->getSize();
  cout << "Number of jets " << n_jets << endl;
  bool JetsInAcceptance[n_jets]; //Mas for keeping track of jets inside acceptance
  bool jetIsBTagged[n_jets]; //Mask for keeping track of b-tagged jets
  int CountingGoodJets=0;
  int CountingBadJets=0;
  float TotalHT=0;

  ScaleFactor jetSF[n_jets];

  TLorentzVector AllJets[n_jets];

  //Loop over jets

  for (int i=0;i<n_jets;++i)
    {
      TLorentzVector *jeti = m_jetMet->getP4(i);
      
      //cout << "Pt of jet " << i << " is " << jeti->Pt() << endl;

      TotalHT+=fabs(jeti->Pt());
      if (!SingleTprime_analysis::isJetSel(jeti)) continue; // apply the pt cut
      AllJets[i].SetPxPyPzE(jeti->Px(),jeti->Py(),jeti->Pz(),jeti->E());      
      if (SingleTprime_analysis::isJetAccepSel(jeti)) {JetsInAcceptance[i]=true; CountingGoodJets++;}
      else {JetsInAcceptance[i]=false;}
      if (SingleTprime_analysis::isJetForwSel(jeti)) {if (!JetsInAcceptance[i]) CountingBadJets++;}

      if (m_isMC) jetSF[i] = m_jetMet->getScaleFactor(i);

      //Finding B-tagged jets with loose requirement based on CSV algorithm
      if ((m_jetMet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVL) 
	{
	  ++m_NBtaggedJets_CSVL;
	  jetIsBTagged[i] = true;
	} 
      else 
	{
	  jetIsBTagged[i] = false;
	}
    }

  // First check the number of jets in and out the acceptance
  cout << CountingGoodJets << " " << CountingBadJets << endl;
  if (CountingGoodJets>=NumberOfGoodJets && CountingBadJets>=NumberOfBadJets) cout << "Good Event" << endl;
  else return 0; 
  
  //cout << "The HT of the event is " << TotalHT << endl;
  m_THT = TotalHT;
  if (n_jets>0) m_jet1pt = AllJets[0].Pt();
  if (n_jets>1) m_jet2pt = AllJets[1].Pt();
  if (n_jets>2) m_jet3pt = AllJets[2].Pt();
  if (n_jets>3) m_jet4pt = AllJets[3].Pt();
  if (n_jets>4) m_jet5pt = AllJets[4].Pt();
  if (n_jets>5) m_jet6pt = AllJets[5].Pt();

  ///////////////////////////////////////////////
  //Reconstructing the Higgs from b-tagged jets//
  ///////////////////////////////////////////////
  int IndexHiggsJets[2]={0,0};
  bool EventWithHiggs=false;
  int HiggsCounter=0;
  float DiffWithHiggMass=0;
  for (int i=0;i<n_jets;++i)
    {
      if (!jetIsBTagged[i] || !JetsInAcceptance[i]) continue;
      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E()); //m_jetMet->getP4(i);
      for (int j=i+1;j<n_jets;++j)
	{
	  if (!jetIsBTagged[j] || !JetsInAcceptance[i]) continue;
	  //TLorentzVector *jeti = m_jetMet->getP4(i);
	  TLorentzVector jetj; jetj.SetPxPyPzE(AllJets[j].Px(),AllJets[j].Py(),AllJets[j].Pz(),AllJets[j].E()); //m_jetMet->getP4(j);
	  TLorentzVector DiBjet = jeti+jetj; //DiBjet->SetPxPyPzE(0,0,0,0);
	  //DiBjet->SetPxPyPzE(jeti->Px()+jetj->Px(),jeti->Py()+jetj->Py(),jeti->Pz()+jetj->Pz(),jeti->E()+jetj->E());
	  if (fabs(DiBjet.M()-HiggsMass)>HiggsMassWindow) continue;
	  //cout << "Mass of jet i: " << jeti.M() << " and j: " << jetj.M() << endl;
	  //cout << "Mass of the b-tagged jets couple " << i << j << " is " << DiBjet.M() << endl;
	  EventWithHiggs=true;
	  if (HiggsCounter==0)
	    {
	      ++HiggsCounter; DiffWithHiggMass=fabs(DiBjet.M()-HiggsMass);
	      IndexHiggsJets[0]=i; IndexHiggsJets[1]=j; ReconstructedHiggs->SetPxPyPzE(jeti.Px()+jetj.Px(),jeti.Py()+jetj.Py(),jeti.Pz()+jetj.Pz(),jeti.E()+jetj.E());
	      FirstHiggsJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E()); SecondHiggsJet->SetPxPyPzE(jetj.Px(),jetj.Py(),jetj.Pz(),jetj.E());
	    }
	  else
	    {
	      ++HiggsCounter;
	      if (DiffWithHiggMass>fabs(DiBjet.M()-HiggsMass)) 
		{
		  IndexHiggsJets[0]=i; IndexHiggsJets[1]=j; 
		  DiffWithHiggMass=fabs(DiBjet.M()-HiggsMass); ReconstructedHiggs->SetPxPyPzE(jeti.Px()+jetj.Px(),jeti.Py()+jetj.Py(),jeti.Pz()+jetj.Pz(),jeti.E()+jetj.E()); 
		  FirstHiggsJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E()); SecondHiggsJet->SetPxPyPzE(jetj.Px(),jetj.Py(),jetj.Pz(),jetj.E());
		}
	    }
	}
    }

  //cout << "Mass of selected Higgs is " << ReconstructedHiggs->M() << " from couple " << IndexHiggsJets[0] << IndexHiggsJets[1] << endl;
  //TLorentzVector BJetCouple;  BJetCouple.SetPxPyPzE(FirstHiggsJet->Px()+SecondHiggsJet->Px(),FirstHiggsJet->Py()+SecondHiggsJet->Py(),FirstHiggsJet->Pz()+SecondHiggsJet->Pz(),FirstHiggsJet->E()+SecondHiggsJet->E());
  //cout << "Mass of selected couple of b jets is " << BJetCouple.M() << endl;

  if (!EventWithHiggs) return 0;
  TLorentzVector FHJ; FHJ.SetPxPyPzE(FirstHiggsJet->Px(), FirstHiggsJet->Py(), FirstHiggsJet->Pz(), FirstHiggsJet->E());
  TLorentzVector SHJ; SHJ.SetPxPyPzE(SecondHiggsJet->Px(), SecondHiggsJet->Py(), SecondHiggsJet->Pz(), SecondHiggsJet->E());
  m_DRHiggsJets=FHJ.DeltaR(SHJ);

  //////////////////////////////////////////////////////////////////
  //Reconstruction of W from full jets collection minus Higgs jets//
  //////////////////////////////////////////////////////////////////
  int IndexWJets[2]={0,0};
  bool EventWithW=false;
  int WCounter=0;
  float DiffWithWMass=0;
  for (int i=0;i<n_jets;++i)
    {
      if (IndexHiggsJets[0]==i || IndexHiggsJets[1]==i || !JetsInAcceptance[i]) continue;
      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
      for (int j=i+1;j<n_jets;++j)
	{
	  if (IndexHiggsJets[0]==j || IndexHiggsJets[1]==j || !JetsInAcceptance[i]) continue;
	  TLorentzVector jetj; jetj.SetPxPyPzE(AllJets[j].Px(),AllJets[j].Py(),AllJets[j].Pz(),AllJets[j].E());
	  TLorentzVector Dijet = jeti+jetj;
	  if (fabs(Dijet.M()-WMass)>WMassWindow) continue;
	  //cout << "Mass of jet i: " << jeti.M() << " and j: " << jetj.M() << endl;
	  //cout << "Mass of the b-tagged jets couple " << i << j << " is " << DiBjet.M() << endl;
	  EventWithW=true;
	  if (WCounter==0)
	    {
	      ++WCounter; DiffWithWMass=fabs(Dijet.M()-WMass);
	      IndexWJets[0]=i; IndexWJets[1]=j; ReconstructedW->SetPxPyPzE(jeti.Px()+jetj.Px(),jeti.Py()+jetj.Py(),jeti.Pz()+jetj.Pz(),jeti.E()+jetj.E());
	      FirstWJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E()); SecondWJet->SetPxPyPzE(jetj.Px(),jetj.Py(),jetj.Pz(),jetj.E());
	    }
	  else
	    {
	      ++WCounter;
	      if (DiffWithWMass>fabs(Dijet.M()-WMass)) 
		{
		  IndexWJets[0]=i; IndexWJets[1]=j; 
		  DiffWithWMass=fabs(Dijet.M()-WMass); ReconstructedW->SetPxPyPzE(jeti.Px()+jetj.Px(),jeti.Py()+jetj.Py(),jeti.Pz()+jetj.Pz(),jeti.E()+jetj.E()); 
		  FirstWJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E()); SecondWJet->SetPxPyPzE(jetj.Px(),jetj.Py(),jetj.Pz(),jetj.E());
		}
	    }
	}
    }

  if (!EventWithW) return 0;
  //cout << "Higgs Jets are: " << IndexHiggsJets[0] << IndexHiggsJets[1] << " W Jets are: " << IndexWJets[0] << IndexWJets[1] << endl;
  TLorentzVector FWJ; FWJ.SetPxPyPzE(FirstWJet->Px(), FirstWJet->Py(), FirstWJet->Pz(), FirstWJet->E());
  TLorentzVector SWJ; SWJ.SetPxPyPzE(SecondWJet->Px(), SecondWJet->Py(), SecondWJet->Pz(), SecondWJet->E());
  m_DRWJets=FHJ.DeltaR(SWJ);  

  /////////////////////////
  //Reconstruction of Top//
  /////////////////////////
  int IndexTopJet=0;
  bool EventWithTop=false;
  int TopCounter=0;
  float DiffWithTopMass=0;
  for (int i=0;i<n_jets;++i)
    {
      if (IndexHiggsJets[0]==i || IndexHiggsJets[1]==i || IndexWJets[0]==i || IndexWJets[1]==i || !JetsInAcceptance[i]) continue;
      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
      TLorentzVector wjet; wjet.SetPxPyPzE(ReconstructedW->Px(),ReconstructedW->Py(),ReconstructedW->Pz(),ReconstructedW->E());
      TLorentzVector Trijet = jeti+wjet;
      if (fabs(Trijet.M()-TopMass)>TopMassWindow) continue;
      //cout << "Mass of jet i: " << jeti.M() << " and j: " << jetj.M() << endl;
      //cout << "Mass of the b-tagged jets couple " << i << j << " is " << DiBjet.M() << endl;
      EventWithTop=true;
      if (TopCounter==0)
	{
	  ++TopCounter; DiffWithTopMass=fabs(Trijet.M()-TopMass);
	  IndexTopJet=i; ReconstructedTop->SetPxPyPzE(jeti.Px()+wjet.Px(),jeti.Py()+wjet.Py(),jeti.Pz()+wjet.Pz(),jeti.E()+wjet.E());
	  TopJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E());
	}
      else
	{
	  ++TopCounter;
	  if (DiffWithTopMass>fabs(Trijet.M()-TopMass)) 
	    {
		  IndexTopJet=i;
		  DiffWithTopMass=fabs(Trijet.M()-TopMass); ReconstructedTop->SetPxPyPzE(jeti.Px()+wjet.Px(),jeti.Py()+wjet.Py(),jeti.Pz()+wjet.Pz(),jeti.E()+wjet.E()); 
		  TopJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E());
	    }
	}
    }

  if (!EventWithTop) return 0;
  //cout << "Higgs Jets are: " << IndexHiggsJets[0] << IndexHiggsJets[1] << " W Jets are: " << IndexWJets[0] << IndexWJets[1] << " Top jet is: " << IndexTopJet << endl;
  TLorentzVector HJ; HJ.SetPxPyPzE(ReconstructedHiggs->Px(), ReconstructedHiggs->Py(), ReconstructedHiggs->Pz(), ReconstructedHiggs->E());
  TLorentzVector TJ; TJ.SetPxPyPzE(ReconstructedTop->Px(), ReconstructedTop->Py(), ReconstructedTop->Pz(), ReconstructedTop->E());
  TLorentzVector WJ; WJ.SetPxPyPzE(ReconstructedW->Px(), ReconstructedW->Py(), ReconstructedW->Pz(), ReconstructedW->E());
  m_DRTopHiggs=HJ.DeltaR(TJ);
  m_DRWHiggs=HJ.DeltaR(WJ);
  m_RelTHT=(ReconstructedHiggs->Pt()+ReconstructedTop->Pt())/TotalHT; //Relative total hadronic energy

  ReconstructedTprime->SetPxPyPzE(ReconstructedHiggs->Px()+ReconstructedTop->Px(),ReconstructedHiggs->Py()+ReconstructedTop->Py(),ReconstructedHiggs->Pz()+ReconstructedTop->Pz(),ReconstructedHiggs->E()+ReconstructedTop->E());

  ///////////////////////////
  //Comparing with MC truth//
  ///////////////////////////

  //cout << "True Tprime components " << TrueTprime->Px() << " " << TrueTprime->Py() << " " << TrueTprime->Pz() << " " << TrueTprime->E() << endl;
  //cout << "Reconstructed Tprime components "<< ReconstructedTprime->Px() << " " << ReconstructedTprime->Py() << " " << ReconstructedTprime->Pz() << " " << ReconstructedTprime->E() << endl;

  if (TrueTprime->Pt()!=0 && ReconstructedTprime->Pt()!=0)
    {
      TLorentzVector TTp; TTp.SetPxPyPzE(TrueTprime->Px(), TrueTprime->Py(), TrueTprime->Pz(), TrueTprime->E());
      TLorentzVector RTp; RTp.SetPxPyPzE(ReconstructedTprime->Px(), ReconstructedTprime->Py(), ReconstructedTprime->Pz(), ReconstructedTprime->E());
      //cout << "DeltaR between true and reconstructed Tprime " << TTp.DeltaR(RTp) << endl;
      if (RTp.DeltaR(TTp)<=m_DRMatching) CorrectTprime=1;
    }
  if (TrueHiggs->Pt()!=0 && ReconstructedHiggs->Pt()!=0)
    {
      TLorentzVector TH; TH.SetPxPyPzE(TrueHiggs->Px(), TrueHiggs->Py(), TrueHiggs->Pz(), TrueHiggs->E());
      TLorentzVector RH; RH.SetPxPyPzE(ReconstructedHiggs->Px(), ReconstructedHiggs->Py(), ReconstructedHiggs->Pz(), ReconstructedHiggs->E());
      if (RH.DeltaR(TH)<=m_DRMatching) CorrectH=1;
    }
  if (FirstTrueHiggsJet->Pt()!=0 && FirstHiggsJet->Pt()!=0)
    {
      TLorentzVector TrHFJ; TrHFJ.SetPxPyPzE(FirstTrueHiggsJet->Px(), FirstTrueHiggsJet->Py(), FirstTrueHiggsJet->Pz(), FirstTrueHiggsJet->E());
      TLorentzVector HFJet; HFJet.SetPxPyPzE(FirstHiggsJet->Px(), FirstHiggsJet->Py(), FirstHiggsJet->Pz(), FirstHiggsJet->E());
      if (HFJet.DeltaR(TrHFJ)<=m_DRMatching) CorrectHiggsJet=1;
    }
  if (SecondTrueHiggsJet->Pt()!=0 && SecondHiggsJet->Pt()!=0)
    {
      TLorentzVector TrHSJ; TrHSJ.SetPxPyPzE(SecondTrueHiggsJet->Px(), SecondTrueHiggsJet->Py(), SecondTrueHiggsJet->Pz(), SecondTrueHiggsJet->E());
      TLorentzVector HSJet; HSJet.SetPxPyPzE(SecondHiggsJet->Px(), SecondHiggsJet->Py(), SecondHiggsJet->Pz(), SecondHiggsJet->E());
      if (HSJet.DeltaR(TrHSJ)<=m_DRMatching) CorrectHiggsJet=2;
    }
  if (TrueW->Pt()!=0 && ReconstructedW->Pt()!=0)
    {
      TLorentzVector TrW; TrW.SetPxPyPzE(TrueW->Px(), TrueW->Py(), TrueW->Pz(), TrueW->E());
      TLorentzVector ReW; ReW.SetPxPyPzE(ReconstructedW->Px(), ReconstructedW->Py(), ReconstructedW->Pz(), ReconstructedW->E());
      if (ReW.DeltaR(TrW)<=m_DRMatching) CorrectW=1;
    }
  if (FirstTrueWJet->Pt()!=0 && FirstWJet->Pt()!=0)
    {
      TLorentzVector TWFJ; TWFJ.SetPxPyPzE(FirstTrueWJet->Px(), FirstTrueWJet->Py(), FirstTrueWJet->Pz(), FirstTrueWJet->E());
      TLorentzVector WFJ; WFJ.SetPxPyPzE(FirstWJet->Px(), FirstWJet->Py(), FirstWJet->Pz(), FirstWJet->E());
      if (WFJ.DeltaR(TWFJ)<=m_DRMatching) CorrectWJet=1;
    }
  if (SecondTrueWJet->Pt()!=0 && SecondWJet->Pt()!=0)
    {
      TLorentzVector TWSJ; TWSJ.SetPxPyPzE(SecondTrueWJet->Px(), SecondTrueWJet->Py(), SecondTrueWJet->Pz(), SecondTrueWJet->E());
      TLorentzVector WSJ; WSJ.SetPxPyPzE(SecondWJet->Px(), SecondWJet->Py(), SecondWJet->Pz(), SecondWJet->E());
      if (WSJ.DeltaR(TWSJ)<=m_DRMatching) CorrectWJet=2;
    }
  if (TrueTop->Pt()!=0 && ReconstructedTop->Pt()!=0)
    {
      TLorentzVector TT; TT.SetPxPyPzE(TrueTop->Px(), TrueTop->Py(), TrueTop->Pz(), TrueTop->E());
      TLorentzVector RT; RT.SetPxPyPzE(ReconstructedTop->Px(), ReconstructedTop->Py(), ReconstructedTop->Pz(), ReconstructedTop->E());
      if (RT.DeltaR(TT)<=m_DRMatching) CorrectTop=1;
    }
  if (TopJet->Pt()!=0 && TopTrueJet->Pt()!=0)
    {
      TLorentzVector TTJ; TTJ.SetPxPyPzE(TopTrueJet->Px(), TopTrueJet->Py(), TopTrueJet->Pz(), TopTrueJet->E());
      TLorentzVector RTJ; RTJ.SetPxPyPzE(TopJet->Px(), TopJet->Py(), TopJet->Pz(), TopJet->E());
      if (RTJ.DeltaR(TTJ)<=m_DRMatching) CorrectTopJet=1;
    }

  // Finally fill the analysis tree 
  //SingleTprime_analysis::fillTree();

  return 0;
}

bool SingleTprime_analysis::isJetSel(TLorentzVector *jet) //Function to select Jets of interest
{  
  double jetpt = jet->Pt();
  double jeteta = jet->Eta();

  if (jetpt<m_jet_Ptcut) { if (fabs(jeteta)>m_jet_EtaMaxcut) return false;}

  return true;
}

bool SingleTprime_analysis::isJetAccepSel(TLorentzVector *jet) //Function to select jets inside the acceptance
{  
  double jeteta = jet->Eta();

  if (fabs(jeteta)<m_jet_EtaAccepcut) return true;

  return false;
}

bool SingleTprime_analysis::isJetForwSel(TLorentzVector *jet) // Function to select jets outside acceptance
{  
  double jeteta = jet->Eta();

  if (fabs(jeteta)>(m_jet_EtaAccepcut-m_jet_OverlapAccep) && fabs(jeteta)<m_jet_EtaMaxcut) return true;

  return false;
}

void SingleTprime_analysis::analyze(const edm::Event& event, const edm::EventSetup& iSetup, PatExtractor& extractor) {
  analyze(iSetup, extractor);
}

void SingleTprime_analysis::analyze(const edm::EventSetup& iSetup, PatExtractor& extractor)
{
  reset(); //Setting all variables to zero

  //Pointer to Extractor objects
  m_jetMet   = std::static_pointer_cast<JetMETExtractor>(extractor.getExtractor("JetMET"));
  m_event    = std::static_pointer_cast<EventExtractor>(extractor.getExtractor("event"));
  if (m_isMC) 
    {
      m_MC = std::static_pointer_cast<MCExtractor>(extractor.getExtractor("MC"));
      MCidentification();
    }

  m_nPU = m_event->nPU();
  m_evt = m_event->n_events();
  evt_num++; cout << evt_num << endl;

  //Performing analysis
  //cout << "Going to selection" << endl;
  SingleTprime_Sel();

  fillTree();
}

  //MC Identification

#define ID_B (5)
#define ID_T (6)
#define ID_H (25)
#define ID_W (24)
#define ID_Tp (6000006)

int SingleTprime_analysis::patIndexToExtractorIndex(int patIndex) const {

  for (int i = 0; i < m_MC->getSize() ; i++) {
    if (m_MC->getPatIndex(i) == patIndex)
      return i;
  }

  return -1;
}

void SingleTprime_analysis::MCidentification()
{

  int n_MC = m_MC->getSize();

  if (!n_MC) return;

  for (int i = 0; i < n_MC ; ++i)
    {

      int motherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(i));
      int grandMotherIndex = -1;
      if (motherIndex != -1) grandMotherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(motherIndex));

      if (abs(m_MC->getType(i)) == ID_H) 
	{
	  TrueHiggs->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}
      
      if (abs(m_MC->getType(i)) == ID_W) 
	{
	  TrueW->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}
      
      if (abs(m_MC->getType(i)) == ID_T) 
	{
	  TrueTop->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (abs(m_MC->getType(i)) == ID_Tp) 
	{
	  TrueTprime->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (m_MC->getType(i) == ID_B && abs(m_MC->getType(motherIndex)) == ID_H) 
	{
	  FirstTrueHiggsJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (m_MC->getType(i) == -1*ID_B && abs(m_MC->getType(motherIndex)) == ID_H) 
	{
	  SecondTrueHiggsJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (m_MC->getType(i) < ID_B && m_MC->getType(i) > 0 && abs(m_MC->getType(motherIndex)) == ID_W) 
	{
	  FirstTrueWJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (m_MC->getType(i) > -1*ID_B && m_MC->getType(i) < 0 && abs(m_MC->getType(motherIndex)) == ID_W) 
	{
	  SecondTrueWJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (abs(m_MC->getType(i)) == ID_B && abs(m_MC->getType(motherIndex)) == ID_T) 
	{
	  TopTrueJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (motherIndex == -1) continue;

      if (true) 
	{
	  std::cout << "Type: " << m_MC->getType(i) << std::endl;
	  std::cout << "Mother type: " << m_MC->getType(motherIndex) << std::endl;
	  if (grandMotherIndex != -1) std::cout << "Grandmother type: " << m_MC->getType(grandMotherIndex) << std::endl;
	}

    }
}

void SingleTprime_analysis::reset()
{
  m_evt = 0;
  m_nPU = 0;

  m_THT = 0.;
  m_jet1pt = 0.;
  m_jet2pt = 0.;
  m_jet3pt = 0.;
  m_jet4pt = 0.;
  m_jet5pt = 0.;
  m_jet6pt = 0.;
  m_NBtaggedJets_CSVL=0;
  m_DRHiggsJets=0.;
  m_DRWJets=0.;
  m_DRTopHiggs=0.;
  m_DRWHiggs=0.;
  m_RelTHT=0.;

  CorrectTprime = 0;
  CorrectH = 0;
  CorrectW = 0;
  CorrectTop = 0;
  CorrectHiggsJet = 0;
  CorrectWJet = 0;
  CorrectTopJet = 0;

  ReconstructedHiggs->SetPxPyPzE(0., 0., 0., 0.);
  ReconstructedW->SetPxPyPzE(0., 0., 0., 0.);
  ReconstructedTop->SetPxPyPzE(0., 0., 0., 0.);
  ReconstructedTprime->SetPxPyPzE(0., 0., 0., 0.);
  TrueHiggs->SetPxPyPzE(0., 0., 0., 0.);
  TrueW->SetPxPyPzE(0., 0., 0., 0.);
  TrueTop->SetPxPyPzE(0., 0., 0., 0.);
  TrueTprime->SetPxPyPzE(0., 0., 0., 0.);
  FirstTrueHiggsJet->SetPxPyPzE(0., 0., 0., 0.);
  SecondTrueHiggsJet->SetPxPyPzE(0., 0., 0., 0.);
  FirstTrueWJet->SetPxPyPzE(0., 0., 0., 0.);
  SecondTrueWJet->SetPxPyPzE(0., 0., 0., 0.);
  TopTrueJet->SetPxPyPzE(0., 0., 0., 0.);
}

// Fill the root tree containing analysis results

void SingleTprime_analysis::fillTree()
{
  m_tree_stp->Fill(); 
}

}

// Register the plugin inside the factory
DEFINE_EDM_PLUGIN(PatExtractorPluginFactory, patextractor::SingleTprime_analysis, "SingleTprime_analysis");
