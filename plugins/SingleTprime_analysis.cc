#include "SingleTprime_analysis.h"

#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

#include "TH2.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

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
int NumberOfBadJets=6;
float THTcut=630;
float HiggsMass=125.0;
float HiggsMassWindow=2000.0;
float WMass=80.3;
float WMassWindow=2000.0;
float TopMass=172.5;
float TopMassWindow=2000.0;
float DeltaRHiggsJets=2.5;
float DeltaRWJets=3.0;

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
  TrueTprimeAcompainingJet = new TLorentzVector(0., 0., 0., 0.);
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
  m_tree_stp->Branch("True_Tprime_Acompaining_Jet", &TrueTprimeAcompainingJet);
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
  m_tree_stp->Branch("DeltaR_TrueHiggs_RecoHiggs",  &m_DRTrueHiggsRecoHiggs   ,"DRTrueHiggsRecoHiggs/F");
  m_tree_stp->Branch("DeltaR_TrueW_RecoW",  &m_DRTrueWRecoW   ,"DRTrueWRecoW/F");
  m_tree_stp->Branch("DeltaR_TrueTop_RecoTop",  &m_DRTrueTopRecoTop   ,"DRTrueTopRecoTop/F");
  m_tree_stp->Branch("DeltaR_TrueTprime_RecoTprime",  &m_DRTrueTprimeRecoTprime   ,"DRTrueTprimeRecoTprime/F");
  m_tree_stp->Branch("DeltaR_TrueFirstHiggsJet_RecoJet",  &m_DRTrueFirstHiggsJetRecoJet   ,"DRTrueFirstHiggsJetRecoJet/F");
  m_tree_stp->Branch("DeltaR_TrueSecondHiggsJet_RecoJet",  &m_DRTrueSecondHiggsJetRecoJet   ,"DRTrueSecondHiggsJetRecoJet/F");
  m_tree_stp->Branch("DeltaR_TrueFirstWJet_RecoJet",  &m_DRTrueFirstWJetRecoJet   ,"DRTrueFirstWJetRecoJet/F");
  m_tree_stp->Branch("DeltaR_TrueSecondWJet_RecoJet",  &m_DRTrueSecondWJetRecoJet   ,"DRTrueSecondWJetRecoJet/F");
  m_tree_stp->Branch("DeltaR_TrueTopJet_RecoJet",  &m_DRTrueTopJetRecoJet   ,"DRTrueTopJetRecoJet/F");
  m_tree_stp->Branch("DeltaR_TrueWJets",  &m_DRTrueWJets   ,"DRTrueWJets/F");
  m_tree_stp->Branch("DeltaR_MatchedWJets",  &m_DRMatchedWJets   ,"DRMatchedWJets/F");
  m_tree_stp->Branch("DeltaPhi_TrueWJets",  &m_DPhiTrueWJets   ,"DRTrueWJets/F");
  m_tree_stp->Branch("DeltaPhi_MatchedWJets",  &m_DPhiMatchedWJets   ,"DRMatchedTopJets/F");
  m_tree_stp->Branch("Number_of_Matched_Higgs_Jets",  &NumbMatchedHiggsJets   ,"NumbMatchedHiggsJets/I");
  m_tree_stp->Branch("Number_of_Matched_W_Jets",  &NumbMatchedWJets   ,"NumbMatchedWJets/I");
  m_tree_stp->Branch("Number_of_Matched_Top_Jets",  &NumbMatchedTopJets   ,"NumbMatchedTopJets/I");
  m_tree_stp->Branch("Sphericity",  &m_Sphericity   ,"Sphericity/F");
  m_tree_stp->Branch("Aplanarity",  &m_Aplanarity   ,"Aplanarity/F");

  // Initialize the analysis parameters using the ParameterSet iConfig
  //int an_option = iConfig.getUntrackedParameter<int>("an_option", 0);
  m_jet_Ptcut = 20;                               // Default val
  m_jet_EtaMaxcut = 4.5;
  m_jet_EtaAccepcut = 2.5;
  m_jet_OverlapAccep = 2.5;
  m_JET_btag_CSVL = 0.244;
  evt_num = 0;
  m_DRMatching=0.3;
  m_DPtMatching=10.0;
}

SingleTprime_analysis::~SingleTprime_analysis(){}

  int SingleTprime_analysis::SingleTprime_Sel() //Main function for the analysis
{
  int n_jets = m_jetMet->getSize();
  //cout << "Number of jets " << n_jets << endl;
  bool JetsInAcceptance[n_jets]; //Mas for keeping track of jets inside acceptance
  bool jetIsBTagged[n_jets]; //Mask for keeping track of b-tagged jets
  int CountingGoodJets=0;
  int CountingBadJets=0;
  float TotalHT=0;

  ScaleFactor jetSF[n_jets];

  TLorentzVector AllJets[n_jets];

  //Loop over jets

  /*int n_MC = m_MC->getSize();

  if (!n_MC) return 0;

  for (int i = 0; i < n_MC ; ++i)
    {
      if (abs(m_MC->getType(i))==0) continue;
      if (abs(m_MC->getType(i)) <= 5 || abs(m_MC->getType(i))==21) 
	{
	  cout << "HERE1" << endl;
	  //if (m_MC->getStatus(i)!=1) continue;
	  TLorentzVector *jetMC= m_MC->getP4(i); 
	  if (jetMC->Pt()==0 || fabs(jetMC->Eta())<=0.01) continue;
	  cout << "HERE2" << endl;
	  cout << jetMC->Pt() << endl;
	  //jeti->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	  if (isJetForwSel(jetMC)) CountingBadJets++;
	  if (isJetAccepSel(jetMC)) CountingGoodJets++;
	}
    }
	
  if (CountingGoodJets>=NumberOfGoodJets && CountingBadJets>=NumberOfBadJets) cout << "Good Event" << endl;
  else return 0; */
	  
  for (int i=0;i<n_jets;++i)
    {
      TLorentzVector *jeti = m_jetMet->getP4(i);
      
      //cout << "Pt of jet " << i << " is " << jeti->Pt() << endl;

      TotalHT+=fabs(jeti->Pt());
      AllJets[i].SetPxPyPzE(jeti->Px(),jeti->Py(),jeti->Pz(),jeti->E());
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
      
      //if (!isJetSel(jeti)) continue; // apply the pt cut
      if (isJetForwSel(jeti)/* && !JetsInAcceptance[i]*/) CountingBadJets++;
      if (isJetAccepSel(jeti)) {JetsInAcceptance[i]=true; CountingGoodJets++;}
      //if (isJetAccepSel(jeti)) cout << "The pt of InJet " << i << " is " << jeti->Pt() << endl;
      else {JetsInAcceptance[i]=false;}      
      //if (isJetForwSel(jeti) && !JetsInAcceptance[i]) cout << "The pt of OutJet " << i << " is " << jeti->Pt() << endl;
    }

  if (CountingGoodJets>=NumberOfGoodJets && CountingBadJets>=NumberOfBadJets) cout << "Good Event" << endl;
  else return 0;

  // First check the number of jets in and out the acceptance
  //cout << CountingGoodJets << " " << CountingBadJets << endl;
  //if (CountingGoodJets>=NumberOfGoodJets && CountingBadJets>=NumberOfBadJets) cout << "Good Event" << endl;
  //else return 0; 
  
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
  //cout << "Entering Higgs reconstruction" << endl;
  int IndexHiggsJets[2]={0,0};
  bool EventWithHiggs=false;
  int HiggsCounter=0;
  float DiffWithHiggMass=0;
  for (int i=0;i<n_jets;++i)
    {
      if (!jetIsBTagged[i] || !JetsInAcceptance[i]) continue;
      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
      for (int j=i+1;j<n_jets;++j)
	{
	  if (!jetIsBTagged[j] || !JetsInAcceptance[i]) continue;
	  TLorentzVector jetj; jetj.SetPxPyPzE(AllJets[j].Px(),AllJets[j].Py(),AllJets[j].Pz(),AllJets[j].E());
	  TLorentzVector DiBjet = jeti+jetj;
	  //DiBjet.SetPxPyPzE(jeti.Px()+jetj.Px(),jeti.Py()+jetj.Py(),jeti.Pz()+jetj.Pz(),jeti.E()+jetj.E());
	  if (fabs(DiBjet.M()-HiggsMass)>HiggsMassWindow) continue;
	  if (jeti.DeltaR(jetj)>DeltaRHiggsJets) continue;
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
  //cout << "Entering W reconstruction" << endl;
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
	  if (jeti.DeltaR(jetj)>DeltaRWJets) continue;
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
  m_DRWJets=FWJ.DeltaR(SWJ);  

  /////////////////////////
  //Reconstruction of Top//
  /////////////////////////
  //cout << "Entering Top reconstruction" << endl;
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

  int MatchedHiggsJets=0;
  int MatchedFirstHiggsJets=0;
  int MatchedSecondHiggsJets=0;
  int MatchedWJets=0;
  int MatchedFirstWJets=0;
  int MatchedSecondWJets=0;
  int MatchedTopJets=0;
  int MatchedJetsIndexes[5]={-1,-1,-1,-1,-1}; //First Higgs Jet, Second Higgs Jet, First W Jet, Second W Jet, Top Jet

  //int n_jets = m_jetMet->getSize();
  for (int i=0;i<n_jets;++i)
    {
      int jetMatch = m_jetMet->getJetMCIndex(i);
      if (jetMatch!=-10) cout << "Jet " << i << " has index of MC particle matched " << jetMatch << " which has pdg " << m_MC->getType(jetMatch) << " mother pdg index " << m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))) << " and grandmother pdg index " << m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))) << endl;
      if (jetMatch!=-10 && fabs(m_MC->getType(jetMatch))==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) ++MatchedHiggsJets;
      if (jetMatch!=-10 && m_MC->getType(jetMatch)==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++MatchedFirstHiggsJets; MatchedJetsIndexes[0]=i;}
      if (jetMatch!=-10 && m_MC->getType(jetMatch)==-5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++MatchedSecondHiggsJets; MatchedJetsIndexes[1]=i;}
      if (jetMatch!=-10 && fabs(m_MC->getType(jetMatch))<=5 && fabs(m_MC->getType(jetMatch))!=0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==24) ++MatchedWJets;
      if (jetMatch!=-10 && fabs(m_MC->getType(jetMatch))<=5 && m_MC->getType(jetMatch)>0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==24) {++MatchedFirstWJets; MatchedJetsIndexes[2]=i;}
      if (jetMatch!=-10 && fabs(m_MC->getType(jetMatch))<=5 && m_MC->getType(jetMatch)<0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==24) {++MatchedSecondWJets; MatchedJetsIndexes[3]=i;}
      if (jetMatch!=-10 && fabs(m_MC->getType(jetMatch))==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==6) {++MatchedTopJets; MatchedJetsIndexes[4]=i;}
      
    }
  NumbMatchedHiggsJets=MatchedHiggsJets;
  NumbMatchedWJets=MatchedWJets;
  NumbMatchedTopJets=MatchedTopJets;

  //Taking into account only events where a the MC indformation was obtained correctly
  if (FirstTrueHiggsJet->Pt()==0 || SecondTrueHiggsJet->Pt()==0 || FirstTrueWJet->Pt()==0 || SecondTrueWJet->Pt()==0 || TopTrueJet->Pt()==0 || TrueW->Pt()==0 || TrueTprimeAcompainingJet->Pt()==0) return 0;

  //Disambiguation with Pt
  TLorentzVector TrHFJ; TrHFJ.SetPxPyPzE(FirstTrueHiggsJet->Px(), FirstTrueHiggsJet->Py(), FirstTrueHiggsJet->Pz(), FirstTrueHiggsJet->E());
  TLorentzVector TrHSJ; TrHSJ.SetPxPyPzE(SecondTrueHiggsJet->Px(), SecondTrueHiggsJet->Py(), SecondTrueHiggsJet->Pz(), SecondTrueHiggsJet->E());
  TLorentzVector TWFJ; TWFJ.SetPxPyPzE(FirstTrueWJet->Px(), FirstTrueWJet->Py(), FirstTrueWJet->Pz(), FirstTrueWJet->E());
  TLorentzVector TWSJ; TWSJ.SetPxPyPzE(SecondTrueWJet->Px(), SecondTrueWJet->Py(), SecondTrueWJet->Pz(), SecondTrueWJet->E());
  TLorentzVector TTJ; TTJ.SetPxPyPzE(TopTrueJet->Px(), TopTrueJet->Py(), TopTrueJet->Pz(), TopTrueJet->E());

  m_DRTrueWJets=TWFJ.DeltaR(TWSJ);
  m_DPhiTrueWJets=fabs(TWFJ.Phi()-TWSJ.Phi());

  if (MatchedFirstHiggsJets>2)
    {
      float DeltaPtMatchedAndTrueHiggsJets=0;
      int CounterDisam=0;
      for (int i=0;i<n_jets;++i)
	{
	  int jetMatch = m_jetMet->getJetMCIndex(i);
	  TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	  if (CounterDisam==0 && jetMatch!=-10 && m_MC->getType(jetMatch)==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++CounterDisam; DeltaPtMatchedAndTrueHiggsJets=fabs(jeti.Pt()-TrHFJ.Pt()); MatchedJetsIndexes[0]=i;}
	  if (CounterDisam!=0 && jetMatch!=-10 && m_MC->getType(jetMatch)==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25 && DeltaPtMatchedAndTrueHiggsJets>fabs(jeti.Pt()-TrHFJ.Pt())) {DeltaPtMatchedAndTrueHiggsJets=fabs(jeti.Pt()-TrHFJ.Pt()); MatchedJetsIndexes[0]=i;}
	}
    }

  if (MatchedSecondHiggsJets>2)
    {
      float DeltaPtMatchedAndTrueHiggsJets=0;
      int CounterDisam=0;
      for (int i=0;i<n_jets;++i)
	{
	  int jetMatch = m_jetMet->getJetMCIndex(i);
	  TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	  if (CounterDisam==0 && jetMatch!=-10 && m_MC->getType(jetMatch)==-5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++CounterDisam; DeltaPtMatchedAndTrueHiggsJets=fabs(jeti.Pt()-TrHSJ.Pt()); MatchedJetsIndexes[1]=i;}
	  if (CounterDisam!=0 && jetMatch!=-10 && m_MC->getType(jetMatch)==-5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25 && DeltaPtMatchedAndTrueHiggsJets>fabs(jeti.Pt()-TrHSJ.Pt())) {DeltaPtMatchedAndTrueHiggsJets=fabs(jeti.Pt()-TrHSJ.Pt()); MatchedJetsIndexes[1]=i;}
	}
    }

  if (MatchedFirstWJets>2)
    {
      float DeltaPtMatchedAndTrueWJets=0;
      int CounterDisam=0;
      for (int i=0;i<n_jets;++i)
	{
	  int jetMatch = m_jetMet->getJetMCIndex(i);
	  TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	  if (CounterDisam==0 && jetMatch!=-10 && m_MC->getType(jetMatch)<=5 && m_MC->getType(jetMatch)>0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++CounterDisam; DeltaPtMatchedAndTrueWJets=fabs(jeti.Pt()-TWFJ.Pt()); MatchedJetsIndexes[2]=i;}
	  if (CounterDisam!=0 && jetMatch!=-10 && m_MC->getType(jetMatch)<=5 && m_MC->getType(jetMatch)>0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25 && DeltaPtMatchedAndTrueWJets>fabs(jeti.Pt()-TWFJ.Pt())) {DeltaPtMatchedAndTrueWJets=fabs(jeti.Pt()-TWFJ.Pt()); MatchedJetsIndexes[2]=i;}
	}
    }

  if (MatchedSecondWJets>2)
    {
      float DeltaPtMatchedAndTrueWJets=0;
      int CounterDisam=0;
      for (int i=0;i<n_jets;++i)
	{
	  int jetMatch = m_jetMet->getJetMCIndex(i);
	  TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	  if (CounterDisam==0 && jetMatch!=-10 && m_MC->getType(jetMatch)>=-5 && m_MC->getType(jetMatch)<0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++CounterDisam; DeltaPtMatchedAndTrueWJets=fabs(jeti.Pt()-TWSJ.Pt()); MatchedJetsIndexes[3]=i;}
	  if (CounterDisam!=0 && jetMatch!=-10 && m_MC->getType(jetMatch)>=-5 && m_MC->getType(jetMatch)<0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25 && DeltaPtMatchedAndTrueWJets>fabs(jeti.Pt()-TWSJ.Pt())) {DeltaPtMatchedAndTrueWJets=fabs(jeti.Pt()-TWSJ.Pt()); MatchedJetsIndexes[3]=i;}
	}
    }

  if (MatchedTopJets>2)
    {
      float DeltaPtMatchedAndTrueTopJets=0;
      int CounterDisam=0;
      for (int i=0;i<n_jets;++i)
	{
	  int jetMatch = m_jetMet->getJetMCIndex(i);
	  TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	  if (CounterDisam==0 && jetMatch!=-10 && fabs(m_MC->getType(jetMatch))==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==6) {++CounterDisam; DeltaPtMatchedAndTrueTopJets=fabs(jeti.Pt()-TTJ.Pt()); MatchedJetsIndexes[4]=i;}
	  if (CounterDisam!=0 && jetMatch!=-10 && fabs(m_MC->getType(jetMatch)==5) && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==6 && DeltaPtMatchedAndTrueTopJets>fabs(jeti.Pt()-TTJ.Pt())) {DeltaPtMatchedAndTrueTopJets=fabs(jeti.Pt()-TTJ.Pt()); MatchedJetsIndexes[4]=i;}
	}
    }

  if (MatchedJetsIndexes[0]!=-1 && MatchedJetsIndexes[1]!=-1 && MatchedJetsIndexes[2]!=-1 && MatchedJetsIndexes[3]!=-1 && MatchedJetsIndexes[4]!=-1)
    {
      m_DRMatchedWJets=AllJets[MatchedJetsIndexes[2]].DeltaR(AllJets[MatchedJetsIndexes[3]]);
      m_DPhiMatchedWJets=fabs(AllJets[MatchedJetsIndexes[2]].Phi()-AllJets[MatchedJetsIndexes[3]].Phi());
      
      int GoodMatchedJets=0;
      bool AllDiffIndex=true;
      
      for (int j=0; j<5; j++) 
	{
	  for (int k=j+1; k<5; k++) {if (MatchedJetsIndexes[j]==MatchedJetsIndexes[k]) AllDiffIndex=false;}
	}
      
      if (AllDiffIndex)
	{
	  for (int j=0; j<5; j++)
	    {
	      if (IndexHiggsJets[0]==MatchedJetsIndexes[j] || IndexHiggsJets[1]==MatchedJetsIndexes[j] || IndexWJets[0]==MatchedJetsIndexes[j] || IndexWJets[1]==MatchedJetsIndexes[j] || IndexTopJet==MatchedJetsIndexes[j]) ++GoodMatchedJets;
	    }
	}
      else cout << "Shared index between decay prodcuts of tprime after disambiguation " << MatchedJetsIndexes[0] << MatchedJetsIndexes[1] << MatchedJetsIndexes[2] << MatchedJetsIndexes[3] << MatchedJetsIndexes[4] << endl;
      
      if (GoodMatchedJets==5) CorrectTprime=1;
      
      int HiggsGoodMatches=0;
      int WGoodMatches=0;
      for (int j=0; j<2; j++)
	{
	  if (IndexHiggsJets[0]==MatchedJetsIndexes[j] || IndexHiggsJets[1]==MatchedJetsIndexes[j]) HiggsGoodMatches++;
	  if (IndexWJets[0]==MatchedJetsIndexes[j+2] || IndexWJets[1]==MatchedJetsIndexes[j+2]) WGoodMatches++;
	}
      
      if (HiggsGoodMatches==2) CorrectH=1;
      if (WGoodMatches==2) CorrectW=1;
      if (IndexTopJet==MatchedJetsIndexes[4]) CorrectTop=1;
    }
  else return 0;

  ///////////////////////////
  //Comparing with MC truth//
  ///////////////////////////
  //cout << "Entering MC truth" << endl;
  
  int JetsIndexes[5]={0,0,0,0,0}; //FirstHiggsJet, SecondHiggsJet, FirstWJet, SecondWJets, TopJet
  float deltaRJets[5]={0,0,0,0,0};
  int Counters[5]={0,0,0,0,0};
  float deltaPTJets[5]={0,0,0,0,0};

  for (int i=0;i<n_jets;++i)
    {
     TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
      if (FirstTrueHiggsJet->Pt()!=0)
	{
	  if (Counters[0]==0) {deltaRJets[0]=jeti.DeltaR(TrHFJ); deltaPTJets[0]=fabs(jeti.Pt()-TrHFJ.Pt()); JetsIndexes[0]=i; ++Counters[0];}
	  if (Counters[0]!=0 && deltaRJets[0]>jeti.DeltaR(TrHFJ)) {deltaRJets[0]=jeti.DeltaR(TrHFJ); deltaPTJets[0]=fabs(jeti.Pt()-TrHFJ.Pt()); JetsIndexes[0]=i;}
	  //cout << "TrHFJ " << jeti.DeltaR(TrHFJ) << endl;
	}
      if (SecondTrueHiggsJet->Pt()!=0)
	{
	  if (Counters[1]==0) {deltaRJets[1]=jeti.DeltaR(TrHSJ); deltaPTJets[1]=fabs(jeti.Pt()-TrHSJ.Pt()); JetsIndexes[1]=i; ++Counters[1];}
	  if (Counters[1]!=0 && deltaRJets[1]>jeti.DeltaR(TrHSJ)) {deltaRJets[1]=jeti.DeltaR(TrHSJ); deltaPTJets[1]=fabs(jeti.Pt()-TrHSJ.Pt()); JetsIndexes[1]=i;}
	  //cout << "TrHSJ " << jeti.DeltaR(TrHSJ) << endl;
	}
      if (FirstTrueWJet->Pt()!=0)
        {
	  if (Counters[2]==0) {deltaRJets[2]=jeti.DeltaR(TWFJ); deltaPTJets[2]=fabs(jeti.Pt()-TWFJ.Pt()); JetsIndexes[2]=i; ++Counters[2];}
	  if (Counters[2]!=0 && deltaRJets[2]>jeti.DeltaR(TWFJ)) {deltaRJets[2]=jeti.DeltaR(TWFJ); deltaPTJets[2]=fabs(jeti.Pt()-TWFJ.Pt()); JetsIndexes[2]=i;}
          //cout << "TWFJ " << jeti.DeltaR(TWFJ) << endl;
        }
      if (SecondTrueWJet->Pt()!=0)
        {
	  if (Counters[3]==0) {deltaRJets[3]=jeti.DeltaR(TWSJ); deltaPTJets[3]=fabs(jeti.Pt()-TWSJ.Pt()); JetsIndexes[3]=i; ++Counters[3];}
	  if (Counters[3]!=0 && deltaRJets[3]>jeti.DeltaR(TWSJ)) {deltaRJets[3]=jeti.DeltaR(TWSJ); deltaPTJets[3]=fabs(jeti.Pt()-TWSJ.Pt()); JetsIndexes[3]=i;}
          //cout << "TWSJ " << jeti.DeltaR(TWSJ) << endl;
        }
      if (TopTrueJet->Pt()!=0)
        {
	  if (Counters[4]==0) {deltaRJets[4]=jeti.DeltaR(TTJ); deltaPTJets[4]=fabs(jeti.Pt()-TTJ.Pt()); JetsIndexes[4]=i; ++Counters[4];}
	  if (Counters[4]!=0 && deltaRJets[4]>jeti.DeltaR(TTJ)) {deltaRJets[4]=jeti.DeltaR(TTJ); deltaPTJets[4]=fabs(jeti.Pt()-TTJ.Pt()); JetsIndexes[4]=i;}
          //cout << "TTJ " << jeti.DeltaR(TTJ) << endl;
        }
    }
  /*cout << "Final DeltaRFHJ " << deltaRJets[0] << " with index " << JetsIndexes[0] << endl;
  cout << "Final DeltaRSHJ " << deltaRJets[1] << " with index " << JetsIndexes[1] << endl;
  cout << "Final DeltaRFWJ " << deltaRJets[2] << " with index " << JetsIndexes[2] << endl;
  cout << "Final DeltaRSWJ " << deltaRJets[3] << " with index " << JetsIndexes[3] << endl;
  cout << "Final DeltaRTJ " << deltaRJets[4] << " with index " << JetsIndexes[4] << endl;
  m_DRTrueFirstHiggsJetRecoJet=deltaRJets[0];
  m_DRTrueSecondHiggsJetRecoJet=deltaRJets[1];
  m_DRTrueFirstWJetRecoJet=deltaRJets[2];
  m_DRTrueSecondWJetRecoJet=deltaRJets[3];
  m_DRTrueTopJetRecoJet=deltaRJets[4];*/

  for (int k=0; k<4; k++)
    {
      for (int j=k+1; j<5; j++)
	{
	  if (JetsIndexes[k]==JetsIndexes[j] && deltaPTJets[k]<deltaPTJets[j])
	    {
	      //cout << "Entering disambiguation" << endl;
	      deltaRJets[j]=0; JetsIndexes[j]=0; Counters[j]=0;
	      for (int i=0;i<n_jets;++i)
		{
		  TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
		  if (j==1)
		    {
		      if (SecondTrueHiggsJet->Pt()!=0)
			{
			  if (Counters[1]==0) {deltaRJets[1]=jeti.DeltaR(TrHSJ); deltaPTJets[1]=fabs(jeti.Pt()-TrHSJ.Pt()); JetsIndexes[1]=i; ++Counters[1];}
			  if (Counters[1]!=0 && deltaRJets[1]>jeti.DeltaR(TrHSJ)) {deltaRJets[1]=jeti.DeltaR(TrHSJ); deltaPTJets[1]=fabs(jeti.Pt()-TrHSJ.Pt()); JetsIndexes[1]=i;}
			  //cout << "TrHSJ " << jeti.DeltaR(TrHSJ) << endl;
			}
		    }
		  else if (j==2)
		    {
		      if (FirstTrueWJet->Pt()!=0)
			{
			  if (Counters[2]==0) {deltaRJets[2]=jeti.DeltaR(TWFJ); deltaPTJets[2]=fabs(jeti.Pt()-TWFJ.Pt()); JetsIndexes[2]=i; ++Counters[2];}
			  if (Counters[2]!=0 && deltaRJets[2]>jeti.DeltaR(TWFJ)) {deltaRJets[2]=jeti.DeltaR(TWFJ); deltaPTJets[2]=fabs(jeti.Pt()-TWFJ.Pt()); JetsIndexes[2]=i;}
			  //cout << "TWFJ " << jeti.DeltaR(TWFJ) << endl;
			}
		    }
		  else if (j==3)
		    {
		      if (SecondTrueWJet->Pt()!=0)
			{
			  if (Counters[3]==0) {deltaRJets[3]=jeti.DeltaR(TWSJ); deltaPTJets[3]=fabs(jeti.Pt()-TWSJ.Pt()); JetsIndexes[3]=i; ++Counters[3];}
			  if (Counters[3]!=0 && deltaRJets[3]>jeti.DeltaR(TWSJ)) {deltaRJets[3]=jeti.DeltaR(TWSJ); deltaPTJets[3]=fabs(jeti.Pt()-TWSJ.Pt()); JetsIndexes[3]=i;}
			  //cout << "TWSJ " << jeti.DeltaR(TWSJ) << endl;
			}
		    }
		  else if (j==4)
		    {
		      if (TopTrueJet->Pt()!=0)
			{
			  if (Counters[4]==0) {deltaRJets[4]=jeti.DeltaR(TTJ); deltaPTJets[4]=fabs(jeti.Pt()-TTJ.Pt()); JetsIndexes[4]=i; ++Counters[4];}
			  if (Counters[4]!=0 && deltaRJets[4]>jeti.DeltaR(TTJ)) {deltaRJets[4]=jeti.DeltaR(TTJ); deltaPTJets[4]=fabs(jeti.Pt()-TTJ.Pt()); JetsIndexes[4]=i;}
			  //cout << "TTJ " << jeti.DeltaR(TTJ) << endl;
			}
		    }
		}
	    }
	  else if (JetsIndexes[k]==JetsIndexes[j] && deltaPTJets[k]>=deltaPTJets[j])
	    {
	      //cout << "Entering disambiguation" << endl;
	      deltaRJets[k]=0; JetsIndexes[k]=0; Counters[k]=0;
	      for (int i=0;i<n_jets;++i)
		{
		  TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
		  if (k==0)
		    {
		      if (FirstTrueHiggsJet->Pt()!=0)
			{
			  if (Counters[0]==0) {deltaRJets[0]=jeti.DeltaR(TrHFJ); deltaPTJets[0]=fabs(jeti.Pt()-TrHFJ.Pt()); JetsIndexes[0]=i; ++Counters[0];}
			  if (Counters[0]!=0 && deltaRJets[0]>jeti.DeltaR(TrHFJ)) {deltaRJets[0]=jeti.DeltaR(TrHFJ); deltaPTJets[0]=fabs(jeti.Pt()-TrHFJ.Pt()); JetsIndexes[0]=i;}
			  //cout << "TrHFJ " << jeti.DeltaR(TrHFJ) << endl;
			}
		    }
		  else if (k==1)
		    {
		      if (SecondTrueHiggsJet->Pt()!=0)
			{
			  if (Counters[1]==0) {deltaRJets[1]=jeti.DeltaR(TrHSJ); deltaPTJets[1]=fabs(jeti.Pt()-TrHSJ.Pt()); JetsIndexes[1]=i; ++Counters[1];}
			  if (Counters[1]!=0 && deltaRJets[1]>jeti.DeltaR(TrHSJ)) {deltaRJets[1]=jeti.DeltaR(TrHSJ); deltaPTJets[1]=fabs(jeti.Pt()-TrHSJ.Pt()); JetsIndexes[1]=i;}
			  //cout << "TrHSJ " << jeti.DeltaR(TrHSJ) << endl;
			}
		    }
		  else if (k==2)
		    {
		      if (FirstTrueWJet->Pt()!=0)
			{
			  if (Counters[2]==0) {deltaRJets[2]=jeti.DeltaR(TWFJ); deltaPTJets[2]=fabs(jeti.Pt()-TWFJ.Pt()); JetsIndexes[2]=i; ++Counters[2];}
			  if (Counters[2]!=0 && deltaRJets[2]>jeti.DeltaR(TWFJ)) {deltaRJets[2]=jeti.DeltaR(TWFJ); deltaPTJets[2]=fabs(jeti.Pt()-TWFJ.Pt()); JetsIndexes[2]=i;}
			  //cout << "TWFJ " << jeti.DeltaR(TWFJ) << endl;
			}
		    }
		  else if (k==3)
		    {
		      if (SecondTrueWJet->Pt()!=0)
			{
			  if (Counters[3]==0) {deltaRJets[3]=jeti.DeltaR(TWSJ); deltaPTJets[3]=fabs(jeti.Pt()-TWSJ.Pt()); JetsIndexes[3]=i; ++Counters[3];}
			  if (Counters[3]!=0 && deltaRJets[3]>jeti.DeltaR(TWSJ)) {deltaRJets[3]=jeti.DeltaR(TWSJ); deltaPTJets[3]=fabs(jeti.Pt()-TWSJ.Pt()); JetsIndexes[3]=i;}
			  //cout << "TWSJ " << jeti.DeltaR(TWSJ) << endl;
			}
		    }
		}
	    }
	}
    } 
  //cout << "Final DeltaRFHJ " << deltaRJets[0] << " with index " << JetsIndexes[0] << endl;
  //cout << "Final DeltaRSHJ " << deltaRJets[1] << " with index " << JetsIndexes[1] << endl;
  //cout << "Final DeltaRFWJ " << deltaRJets[2] << " with index " << JetsIndexes[2] << endl;
  //cout << "Final DeltaRSWJ " << deltaRJets[3] << " with index " << JetsIndexes[3] << endl;
  //cout << "Final DeltaRTJ " << deltaRJets[4] << " with index " << JetsIndexes[4] << endl;
  m_DRTrueFirstHiggsJetRecoJet=deltaRJets[0];
  m_DRTrueSecondHiggsJetRecoJet=deltaRJets[1];
  m_DRTrueFirstWJetRecoJet=deltaRJets[2];
  m_DRTrueSecondWJetRecoJet=deltaRJets[3];
  m_DRTrueTopJetRecoJet=deltaRJets[4];

  /////////////////////////////////////////////
  //Extracting info for Sphericity/Aplanarity//
  /////////////////////////////////////////////
  TLorentzVector ForDecayRestFrame;
  for (int i=0;i<n_jets;++i)
    {
      TLorentzVector *jeti = m_jetMet->getP4(i);
      TLorentzVector temp;
      if (i==0) ForDecayRestFrame.SetPxPyPzE(jeti->Px(),jeti->Py(),jeti->Pz(),jeti->E());
      else {temp=ForDecayRestFrame; ForDecayRestFrame.SetPxPyPzE(jeti->Px()+temp.Px(),jeti->Py()+temp.Py(),jeti->Pz()+temp.Pz(),jeti->E()+temp.E());}
    }

  double Sxx=0;
  double Sxy=0;
  double Sxz=0;
  double Syy=0;
  double Syz=0;
  double Szz=0;
  double Sxxnum=0;
  double Sxynum=0;
  double Sxznum=0;
  double Syynum=0;
  double Syznum=0;
  double Szznum=0;
  double Sden=0;
    
  for (int i=0;i<n_jets;++i)
    {
      TLorentzVector *jeti = m_jetMet->getP4(i);
      Sxxnum+=jeti->Px()*jeti->Px();
      Sxynum+=jeti->Px()*jeti->Py();
      Sxznum+=jeti->Px()*jeti->Pz();
      Syynum+=jeti->Py()*jeti->Py();
      Syznum+=jeti->Py()*jeti->Pz();
      Szznum+=jeti->Pz()*jeti->Pz();
      Sden+=jeti->P()*jeti->P();
    }
  
  Sxx=Sxxnum/Sden;
  Sxy=Sxynum/Sden;
  Sxz=Sxznum/Sden;
  Syy=Syynum/Sden;
  Syz=Syznum/Sden;
  Szz=Szznum/Sden;
  const int N = 3;
  double s[N*N] = {
    Sxx, Sxy, Sxz,
    Sxy, Syy, Syz,
    Sxz, Syz, Szz
  };
  TMatrixDSym S(N, s);
  TMatrixDSymEigen SEigen(S);
  TVectorD eigenval = SEigen.GetEigenValues();
  
  m_Sphericity=(eigenval[1]+eigenval[2])*3./2.;
  m_Aplanarity=eigenval[2]*3./2.;

  // Finally fill the analysis tree 
  //SingleTprime_analysis::fillTree();

  return 1;
}

bool SingleTprime_analysis::isJetSel(TLorentzVector *jet) //Function to select Jets of interest
{  
  double jetpt = jet->Pt();
  double jeteta = jet->Eta();

  if (jetpt<m_jet_Ptcut || fabs(jeteta)>m_jet_EtaMaxcut) return false;
  
  return true;
}

bool SingleTprime_analysis::isJetAccepSel(TLorentzVector *jet) //Function to select jets inside the acceptance
{  
  double jetpt = jet->Pt();
  double jeteta = jet->Eta();

  if (jetpt>m_jet_Ptcut && fabs(jeteta)<m_jet_EtaAccepcut) return true;

  return false;
}

bool SingleTprime_analysis::isJetForwSel(TLorentzVector *jet) // Function to select jets outside acceptance
{  
  double jetpt = jet->Pt();
  double jeteta = jet->Eta();

  if (jetpt>(m_jet_Ptcut-10) && fabs(jeteta)>=(m_jet_EtaAccepcut-m_jet_OverlapAccep) && fabs(jeteta)<m_jet_EtaMaxcut) return true;

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
  int ToFillTree = SingleTprime_Sel();

  if (ToFillTree==1) fillTree();
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
            
      if (abs(m_MC->getType(i)) == ID_T) 
	{
	  TrueTop->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (abs(m_MC->getType(i)) == ID_Tp) 
	{
	  TrueTprime->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	  for (int j = 0; j < n_MC ; ++j)
	    {
	      
	      int motherIndex1 = patIndexToExtractorIndex(m_MC->getMom1Index(j));
	      int grandMotherIndex1 = -1;
	      if (motherIndex1 != -1) grandMotherIndex1 = patIndexToExtractorIndex(m_MC->getMom1Index(motherIndex1));
	     
	      if (abs(m_MC->getType(j)) < ID_B && abs(m_MC->getType(j))>0 && abs(m_MC->getType(motherIndex)) < ID_B && m_MC->getType(motherIndex)==m_MC->getType(motherIndex1)) 
		{
		  TrueTprimeAcompainingJet->SetPxPyPzE(m_MC->getPx(j),m_MC->getPy(j),m_MC->getPz(j),m_MC->getE(j));
		}
	      
	    }
	}

      if (motherIndex == -1) continue;
      TLorentzVector CurrentMCParticle; CurrentMCParticle.SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
      if (CurrentMCParticle.Pt()==0) continue;

      if (abs(m_MC->getType(i)) == ID_W && abs(m_MC->getType(motherIndex)) == ID_T) 
	{
	  TrueW->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}
      
      if (m_MC->getType(i) == ID_B && abs(m_MC->getType(motherIndex)) == ID_H) 
	{
	  FirstTrueHiggsJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (m_MC->getType(i) == -1*ID_B && abs(m_MC->getType(motherIndex)) == ID_H) 
	{
	  SecondTrueHiggsJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (m_MC->getType(i) <= ID_B && m_MC->getType(i) > 0 && abs(m_MC->getType(motherIndex)) == ID_W) 
	{
	  FirstTrueWJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (m_MC->getType(i) >= -1*ID_B && m_MC->getType(i) < 0 && abs(m_MC->getType(motherIndex)) == ID_W) 
	{
	  SecondTrueWJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (abs(m_MC->getType(i)) == ID_B && abs(m_MC->getType(motherIndex)) == ID_T) 
	{
	  TopTrueJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (false) 
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
  m_Sphericity=0.;
  m_Aplanarity=0.;
  m_DRTrueHiggsRecoHiggs=0.;
  m_DRTrueWRecoW=0.;
  m_DRTrueTopRecoTop=0.;
  m_DRTrueTprimeRecoTprime=0.;
  m_DRTrueFirstHiggsJetRecoJet=0;
  m_DRTrueSecondHiggsJetRecoJet=0;
  m_DRTrueFirstWJetRecoJet=0;
  m_DRTrueSecondWJetRecoJet=0; 
  m_DRTrueTopJetRecoJet=0;

  m_DRTrueWJets=0;   
  m_DRMatchedWJets=0;
  m_DPhiTrueWJets=0;
  m_DPhiMatchedWJets=0;

  CorrectTprime = 0;
  CorrectH = 0;
  CorrectW = 0;
  CorrectTop = 0;
  CorrectHiggsJet = 0;
  CorrectWJet = 0;
  CorrectTopJet = 0;
  NumbMatchedHiggsJets = 0;
  NumbMatchedWJets = 0;
  NumbMatchedTopJets = 0;

  ReconstructedHiggs->SetPxPyPzE(0., 0., 0., 0.);
  ReconstructedW->SetPxPyPzE(0., 0., 0., 0.);
  ReconstructedTop->SetPxPyPzE(0., 0., 0., 0.);
  ReconstructedTprime->SetPxPyPzE(0., 0., 0., 0.);
  TrueHiggs->SetPxPyPzE(0., 0., 0., 0.);
  TrueW->SetPxPyPzE(0., 0., 0., 0.);
  TrueTop->SetPxPyPzE(0., 0., 0., 0.);
  TrueTprime->SetPxPyPzE(0., 0., 0., 0.);
  TrueTprimeAcompainingJet->SetPxPyPzE(0., 0., 0., 0.);
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
