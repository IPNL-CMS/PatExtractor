#include "SingleTprime.h"

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
#include "Extractors/PatExtractor/interface/JetMETExtractor.h"
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

using namespace std;

SingleTprime::SingleTprime(const edm::ParameterSet& iConfig): Plugin(iConfig)
{
  // Initialize the analysis parameters using the ParameterSet iConfig
  int an_option = iConfig.getUntrackedParameter<int>("an_option", 0);
  
  //Histograms
  //Distance = new TH1F("Distance",  "Distance"  , 100, 0, 10);
  Cut0 = new TH1F("Cut0",  "Cut0"  , 150, 0, 1500);
  Cut1 = new TH1F("Cut1",  "Cut1"  , 150, 0, 1500);
  Cut2 = new TH1F("Cut2",  "Cut2"  , 150, 0, 1500);
  Cut3 = new TH1F("Cut3",  "Cut3"  , 150, 0, 1500);
  Cut4 = new TH1F("Cut4",  "Cut4"  , 150, 0, 1500);
  Cut5 = new TH1F("Cut5",  "Cut5"  , 150, 0, 1500);
  Cut6 = new TH1F("Cut6",  "Cut6"  , 150, 0, 1500);
  Cut7 = new TH1F("Cut7",  "Cut7"  , 150, 0, 1500);
  Cut8 = new TH1F("Cut8",  "Cut8"  , 150, 0, 1500);
  Cut9 = new TH1F("Cut9",  "Cut9"  , 150, 0, 1500);
  Cut10 = new TH1F("Cut10",  "Cut10"  , 150, 0, 1500);
  Cut11 = new TH1F("Cut11",  "Cut11"  , 150, 0, 1500);
  Cut12 = new TH1F("Cut12",  "Cut12"  , 150, 0, 1500);
  Cut13 = new TH1F("Cut13",  "Cut13"  , 150, 0, 1500);
  Cut14 = new TH1F("Cut14",  "Cut14"  , 150, 0, 1500);
    //////////////////////////////////
  //Histograms for final cosmetics//
  //////////////////////////////////
  Jet1_pt       = new TH1F("Jet1_pt",  "Jet1_pt"  , 100, 0, 1000);               //After cut 0
  Jet2_pt       = new TH1F("Jet2_pt",  "Jet2_pt"  , 100, 0, 1000);               //After cut 0
  Jet3_pt       = new TH1F("Jet3_pt",  "Jet3_pt"  , 100, 0, 1000);               //After cut 0
  Jet4_pt       = new TH1F("Jet4_pt",  "Jet4_pt"  , 100, 0, 1000);               //After cut 0
  Jet5_pt       = new TH1F("Jet5_pt",  "Jet5_pt"  , 100, 0, 1000);               //After cut 0
  Jet6_pt       = new TH1F("Jet6_pt",  "Jet6_pt"  , 100, 0, 1000);               //After cut 0
  THT           = new TH1F("THT",  "THT"  , 150, 0, 1500);                       //After cut 2
  DeltaRH       = new TH1F("DeltaRH",  "DeltaRH"  , 100, 0, 8);                  //After cut 4
  MHDR          = new TH2F("MHDR",  "MHDR"  , 200, 0, 200, 100, 0, 8);           //After cut 4
  MAroundHiggs  = new TH1F("MAroundHiggs",  "MAroundHiggs"  , 200, 0, 200);      //After cut 4
  MAroundW      = new TH1F("MAroundW",  "MAroundW"  , 200, 0, 200);              //After cut 4
  MAroundTop    = new TH1F("MAroundTop",  "MAroundTop"  , 250, 0, 250);          //After cut 4  
  HptToppt      = new TH2F("HptToppt",  "HptToppt"  , 500, 0, 500, 500, 0, 500); //After cut 5
  DRWDRH        = new TH1F("DRWDRH", "DRWDRH", 100, 0, 8);                       //After cut 6
  DPHDPT        = new TH2F("DPHDPT",  "DPHDPT"  , 100, 0, 7, 100, 0, 7);       //After cut 7
  jet_num       = new TH1F("jet_num",  "jet_num"  , 20, 0, 20);                  //After cut 8
  DPHDPW        = new TH2F("DPHDPW",  "DPHDPW"  , 100, 0, 7, 100, 0, 7);       //After cut 9
  HiggsJetsMass = new TH1F("HiggsJetsMass", "HiggsJetsMass", 200, 0, 200);       //After cut 10
  RelHT         = new TH1F("RelHT",  "RelHT"  , 20, 0, 1);                       //After cut 11
  AplanarityHad = new TH1F("AplanarityHad",  "AplanarityHad"  , 100, 0, 0.5);    //After cut 12
  TprimeMassVsHiggsMass = new TH2F("TprimeMassVsHiggsMass",  "TprimeMassVsHiggsMass"  , 150, 0, 1500, 200, 0, 200);      
  TprimeMassVsDRTopHiggs = new TH2F("TprimeMassVsDRTopHiggs",  "TprimeMassVsDRTopHiggs"  , 150, 0, 1500, 100, 0, 8);
  TprimeMassVsTprimePT = new TH2F("TprimeMassVsTprimePT",  "TprimeMassVsTprimePT"  , 150, 0, 1500, 1000, 0, 1000);
  
}

int SingleTprime::JetSel()
{
  AllJetsPt       = 0.;
  m_selJetsIds.clear();

  int n_jet = m_jetMet->getSize();

  if (! n_jet)
    return 8;

  ScaleFactor jetSF[n_jet];
  bool jetIsBTagged[n_jet];

  for (int i = 0; i < n_jet; i++)
  {
    TLorentzVector *jetP = m_jetMet->getP4(i);

    if (fabs(jetP->Pt()) < m_JET_Pt_min) continue;
    if (fabs(jetP->Eta()) >  m_JET_Eta_max) continue;

    m_mtt_JetEta[m_mtt_NJets] = jetP->Eta();
    m_mtt_JetPt[m_mtt_NJets]  = jetP->Pt();
    if (m_isMC)
      jetSF[m_mtt_NJets] = m_jetMet->getScaleFactor(i);

    ++m_mtt_NJets;

    AllJetsPt += jetP->Pt();

    if (m_mtt_NJets < 9) // Count the number of btagged jets in the selected jets
    {
      m_selJetsIds.push_back(i);
      if ((m_jetMet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVL)
        ++m_mtt_NBtaggedJets_CSVL;

      if ((m_jetMet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVM) {
        ++m_mtt_NBtaggedJets_CSVM;
        jetIsBTagged[i] = true;
      } else {
        jetIsBTagged[i] = false;
      }

      if ((m_jetMet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVT)
        ++m_mtt_NBtaggedJets_CSVT;
      if ((m_jetMet->getJetBTagProb_TCHP(i)) > m_JET_btag_TCHPT)
        ++m_mtt_NBtaggedJets_TCHPT;
    }

    if (m_mtt_NJets == 1) m_mtt_1stjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
    if (m_mtt_NJets == 2) m_mtt_2ndjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
    if (m_mtt_NJets == 3) m_mtt_3rdjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
    if (m_mtt_NJets == 4) m_mtt_4thjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
  }

  // We need at least 4 good jets
  if (m_mtt_NJets < 4)
    return 8;

  if (m_isMC) {
    // Get scale factors
    double sf = 1;
    double sf_error = 0;

    if (m_mtt_NBtaggedJets_CSVM == 1) {

      ScaleFactor sf_jet = {1, 0, 0};

      for (int i = 0; i < m_mtt_NJets; i++)
      {
        if (jetIsBTagged[i]) {
          sf_jet = jetSF[i];
          break;
        }
      }

      // We have one 1 b-tagged jet and not the other one
      // e = e_btag * (1 - e_btag) * 2
      sf = sf_jet.getValue() * ((1. - m_b_tagging_efficiency * sf_jet.getValue()) / (1. - m_b_tagging_efficiency));
      sf_error = ((1 - 2 * m_b_tagging_efficiency * sf_jet.getValue()) / (1 - m_b_tagging_efficiency)) * sf_jet.getErrorHigh();
      sf_error *= sf_error; // square

    } else if (m_mtt_NBtaggedJets_CSVM >= 2) {

      ScaleFactor sf_jet1 = {1, 0, 0};
      ScaleFactor sf_jet2 = {1, 0, 0};
      bool gotJet1 = false;

      for (int i = 0; i < m_mtt_NJets; i++)
      {
        if (jetIsBTagged[i]) {
          if (! gotJet1) {
            gotJet1 = true;
            sf_jet1 = jetSF[i];
          } else {
            sf_jet2 = jetSF[i];
            break;
          }
        }
      }

      sf = sf_jet1.getValue() * sf_jet2.getValue();
      sf_error = sf_jet2.getValue() * sf_jet2.getValue() * sf_jet1.getErrorHigh() * sf_jet1.getErrorHigh() +
        sf_jet1.getValue() * sf_jet1.getValue() * sf_jet2.getErrorHigh() * sf_jet2.getErrorHigh();
    }

    m_weight *= sf;

    double squared_sf_error = sf_error * sf_error;
    m_weight_error_low += squared_sf_error;
    m_weight_error_high += squared_sf_error;
  }

  return 1;
}

SingleTprime::::analyze(const edm::EventSetup& iSetup, PatExtractor& extractor)
{
  reset();

  //Getting links to extracteor objects
  m_vertex   = std::static_pointer_cast<VertexExtractor>(extractor.getExtractor("vertex"));

  m_jetMet   = std::static_pointer_cast<JetMETExtractor>(extractor.getExtractor("JetMET"));
  
  m_event    = std::static_pointer_cast<EventExtractor>(extractor.getExtractor("event"));

  std::shared_ptr<HLTExtractor> HLT = std::static_pointer_cast<HLTExtractor>(extractor.getExtractor("HLT"));
  m_trigger_passed = HLT->isTriggerFired();

  if (m_isMC)
  {
    m_MC = std::static_pointer_cast<MCExtractor>(extractor.getExtractor("MC"));
    MCidentification();
  }

  m_nPU = m_event->nPU();//Numbre of PileUp in each event

  

}

// Register the plugin inside the factory
DEFINE_EDM_PLUGIN(PatExtractorPluginFactory, SingleTprime,  "SingleTprime");
