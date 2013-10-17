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

using namespace std;

namespace patextractor {

SingleTprime_analysis::SingleTprime_analysis(const edm::ParameterSet& iConfig): Plugin(iConfig)
{
  std::cout << "Entering SingleTprime analysis" << std::endl;

  /// Tree definition containing your analysis results

  m_tree_stp = new TTree("stp","Single Analysis info");  
  
  /// Branches definition

  m_tree_stp->Branch("evt",         &m_evt    ,"evt/I");      // Simple evt number or event ID
  m_tree_stp->Branch("THT", &m_THT,"THT/F");
  m_tree_stp->Branch("jet1_pt",  &m_jet1pt   ,"jet1_pt/F");
  m_tree_stp->Branch("jet2_pt",  &m_jet2pt   ,"jet2_pt/F");
  m_tree_stp->Branch("jet3_pt",  &m_jet3pt   ,"jet3_pt/F");
  m_tree_stp->Branch("jet4_pt",  &m_jet4pt   ,"jet4_pt/F");
  m_tree_stp->Branch("jet5_pt",  &m_jet5pt   ,"jet5_pt/F");
  m_tree_stp->Branch("jet6_pt",  &m_jet6pt   ,"jet6_pt/F");
  m_tree_stp->Branch("Number_CSVLbtagged_jets",  &m_NBtaggedJets_CSVL   ,"m_NBtaggedJets_CSVL/I");

  // Initialize the analysis parameters using the ParameterSet iConfig
  //int an_option = iConfig.getUntrackedParameter<int>("an_option", 0);
  m_jet_Ptcut = 20;                               // Default val
  m_jet_EtaMaxcut = 5.0;
  m_jet_EtaAccepcut = 2.5;
  m_jet_MultInAcceptance = 5;
  m_jet_MultOutAcceptance = 1;
  m_JET_btag_CSVL = 0.244;
  //isMC=true;
}

SingleTprime_analysis::~SingleTprime_analysis(){}

int SingleTprime_analysis::SingleTprime_Sel() 
{
  int n_jets = m_jetMet->getSize();
  cout << "Number of jets " << n_jets << endl;
  bool JetsInAcceptance[n_jets]; //Mas for keeping track of jets inside acceptance
  bool jetIsBTagged[n_jets]; //Mask for keeping track of b-tagged jets
  int CountingGoodJets=0;
  int CountingBadJets=0;
  float TotalHT=0;

  //Leading jets
  TLorentzVector *LeadJets[6]; 

  //Loop over jets

  for (int i=0;i<n_jets;++i)
    {
      TLorentzVector *jeti = m_jetMet->getP4(i);
      if (i<6) LeadJets[i] = jeti;
      
      cout << "Pt of jet " << i << " is " << jeti->Pt() << endl;

      TotalHT+=fabs(jeti->Pt());
      if (!SingleTprime_analysis::isJetSel(jeti)) continue; // apply the pt cut
      if (SingleTprime_analysis::isJetAccepSel(jeti)) {JetsInAcceptance[i]=true; CountingGoodJets++;}
      else {JetsInAcceptance[i]=false;}
      if (SingleTprime_analysis::isJetForwSel(jeti)) CountingBadJets++;

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
  
  //m_evt = evtnum;  
  cout << "The HT of the event is " << TotalHT << endl;
  m_THT = TotalHT;
  if (n_jets>0) m_jet1pt = LeadJets[0]->Pt();
  if (n_jets>1) m_jet2pt = LeadJets[1]->Pt();
  if (n_jets>2) m_jet3pt = LeadJets[2]->Pt();
  if (n_jets>3) m_jet4pt = LeadJets[3]->Pt();
  if (n_jets>4) m_jet5pt = LeadJets[4]->Pt();
  if (n_jets>5) m_jet6pt = LeadJets[5]->Pt();

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

  if (fabs(jeteta)>m_jet_EtaAccepcut && fabs(jeteta)<m_jet_EtaMaxcut) return true;

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
  if (m_isMC) m_MC = std::static_pointer_cast<MCExtractor>(extractor.getExtractor("MC"));

  //Performing analysis
  //cout << "Going to selection" << endl;
  SingleTprime_Sel();

  fillTree();
}

void SingleTprime_analysis::reset()
{
  m_evt = 0;

  m_THT = 0.;
  m_jet1pt = 0.;
  m_jet2pt = 0.;
  m_jet3pt = 0.;
  m_jet4pt = 0.;
  m_jet5pt = 0.;
  m_jet6pt = 0.;
  m_NBtaggedJets_CSVL=0;
}

// Fill the root tree containing analysis results

void SingleTprime_analysis::fillTree()
{
  m_tree_stp->Fill(); 
}

}

// Register the plugin inside the factory
DEFINE_EDM_PLUGIN(PatExtractorPluginFactory, patextractor::SingleTprime_analysis, "SingleTprime_analysis");
