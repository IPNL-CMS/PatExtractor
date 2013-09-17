#include <Extractors/PatExtractor/plugins/mtt_analysis.h>

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

using namespace std;

namespace patextractor {

mtt_analysis::mtt_analysis(const edm::ParameterSet& cmsswSettings):
  Plugin(cmsswSettings),
  maxNrIter_                        (cmsswSettings.getParameter<unsigned>     ("maxNrIter"           )),
  maxDeltaS_                        (cmsswSettings.getParameter<double>       ("maxDeltaS"           )),
  maxF_                             (cmsswSettings.getParameter<double>       ("maxF"                )),
  jetParam_                         (cmsswSettings.getParameter<unsigned>     ("jetParametrisation"  )),
  lepParam_                         (cmsswSettings.getParameter<unsigned>     ("lepParametrisation"  )),
  metParam_                         (cmsswSettings.getParameter<unsigned>     ("metParametrisation"  )),
  constraints_                      (cmsswSettings.getParameter<std::vector<unsigned> >("constraints")),
  mW_                               (cmsswSettings.getParameter<double>       ("mW"                  )),
  mTop_                             (cmsswSettings.getParameter<double>       ("mTop"                )),
  jetEnergyResolutionScaleFactors_  (cmsswSettings.getParameter<std::vector<double> >("jetEnergyResolutionScaleFactors")),
  jetEnergyResolutionEtaBinning_    (cmsswSettings.getParameter<std::vector<double> >("jetEnergyResolutionEtaBinning"))
{
 
  m_lepTopP4_AfterChi2 = new TLorentzVector(0., 0., 0., 0.);
  m_hadTopP4_AfterChi2 = new TLorentzVector(0., 0., 0., 0.);

  reset();

  /// Tree definition
  m_tree_Mtt = new TTree("Mtt", "Analysis info");

  /// Branches definition

  m_tree_Mtt->Branch("MC_channel"         , &m_MC_channel             , "MC_channel/I");
  m_tree_Mtt->Branch("MC_mtt"             , &m_MC_mtt                 , "MC_mtt/F");
  m_tree_Mtt->Branch("MC_nPU"             , &m_nPU                    , "m_nPU/I");

  // Indexes of gen particles inside the MC collection. Only valid for semi-lept events
  m_tree_Mtt->Branch("MC_leptonIndex"     , &m_leptonIndex            , "MC_leptonIndex/I");
  m_tree_Mtt->Branch("MC_neutrinoIndex"   , &m_neutrinoIndex          , "MC_neutrinoIndex/I");
  m_tree_Mtt->Branch("MC_leptonicTopIndex", &m_leptonicTopIndex       , "MC_leptonicTopIndex/I");
  m_tree_Mtt->Branch("MC_leptonicBIndex"  , &m_leptonicBIndex         , "MC_leptonicBIndex/I");

  m_tree_Mtt->Branch("MC_hadronicBIndex"  , &m_hadronicBIndex         , "MC_hadronicBIndex/I");
  m_tree_Mtt->Branch("MC_hadronicFirstJetIndex" , &m_firstJetIndex    , "MC_hadronicFirstJetIndex/I");
  m_tree_Mtt->Branch("MC_hadronicSecondJetIndex", &m_secondJetIndex   , "MC_hadronicSecondJetIndex/I");

  m_tree_Mtt->Branch("MC_hadronicWMass"   , &m_MC_hadronicWMass       , "MC_hadronicWMass/F");
  m_tree_Mtt->Branch("MC_leptonicWMass"   , &m_MC_leptonicWMass       , "MC_leptonicWMass/F");
  m_tree_Mtt->Branch("MC_hadronicTopMass" , &m_MC_hadronicTopMass     , "MC_hadronicTopMass/F");
  m_tree_Mtt->Branch("MC_leptonicTopMass" , &m_MC_leptonicTopMass     , "MC_leptonicTopMass/F");
  m_tree_Mtt->Branch("MC_pt_tt"           , &m_MC_pt_tt               , "MC_pt_tt/F");
  m_tree_Mtt->Branch("MC_eta_tt"          , &m_MC_eta_tt              , "MC_eta_tt/F");
  m_tree_Mtt->Branch("MC_beta_tt"         , &m_MC_beta_tt             , "MC_beta_tt/F");

  m_tree_Mtt->Branch("nGoodMuons"         , &m_mtt_NGoodMuons         , "nGoodMuons/I");
  m_tree_Mtt->Branch("nLooseGoodMuons"    , &m_mtt_NLooseGoodMuons    , "nLooseGoodMuons/I");
  m_tree_Mtt->Branch("muonPt"             , &m_mtt_MuonPt             , "muonPt[nGoodMuons]/F");
  m_tree_Mtt->Branch("muonEta"            , &m_mtt_MuonEta            , "muonEta[nGoodMuons]/F");
  m_tree_Mtt->Branch("2DDrMin"            , &m_mtt_2DDrMin            , "2DDrMin[nGoodMuons]/F");
  m_tree_Mtt->Branch("2DpTrel"            , &m_mtt_2DpTrel            , "2DpTrel[nGoodMuons]/F");
  m_tree_Mtt->Branch("muRelIso"           , &m_mtt_MuRelIso           , "muRelIso[nGoodMuons]/F");

  m_tree_Mtt->Branch("nGoodElectrons"     , &m_mtt_NGoodElectrons     , "nGoodElectrons/I");
  m_tree_Mtt->Branch("electronPt"         , &m_mtt_ElectronPt         , "electronPt[nGoodElectrons]/F");
  m_tree_Mtt->Branch("electronEta"        , &m_mtt_ElectronEta        , "electronEta[nGoodElectrons]/F");
  m_tree_Mtt->Branch("elRelIso"           , &m_mtt_ElRelIso           , "elRelIso[nGoodElectrons]/F");
  m_tree_Mtt->Branch("hyperTight1MC"      , &m_mtt_HyperTight1MC      , "hyperTight1MC[nGoodElectrons]/I");


  m_tree_Mtt->Branch("1stjetpt"           , &m_mtt_1stjetpt           , "1stjetpt/F");
  m_tree_Mtt->Branch("2ndjetpt"           , &m_mtt_2ndjetpt           , "2ndjetpt/F");
  m_tree_Mtt->Branch("3rdjetpt"           , &m_mtt_3rdjetpt           , "3rdjetpt/F");
  m_tree_Mtt->Branch("4thjetpt"           , &m_mtt_4thjetpt           , "4thjetpt/F");

  m_tree_Mtt->Branch("nJets"              , &m_mtt_NJets              , "nJets/I");
  m_tree_Mtt->Branch("jetEta"             , &m_mtt_JetEta             , "jetEta[nJets]/F");
  m_tree_Mtt->Branch("jetPt"              , &m_mtt_JetPt              , "jetPt[nJets]/F");

  m_tree_Mtt->Branch("nBtaggedJets_TCHPT" , &m_mtt_NBtaggedJets_TCHPT    , "nBtaggedJets_TCHPT/I");
  m_tree_Mtt->Branch("nBtaggedJets_CSVL" ,  &m_mtt_NBtaggedJets_CSVL    , "nBtaggedJets_CSVL/I");
  m_tree_Mtt->Branch("nBtaggedJets_CSVM" ,  &m_mtt_NBtaggedJets_CSVM    , "nBtaggedJets_CSVM/I");
  m_tree_Mtt->Branch("nBtaggedJets_CSVT" ,  &m_mtt_NBtaggedJets_CSVT    , "nBtaggedJets_CSVT/I");

  m_tree_Mtt->Branch("MET"                , &m_mtt_MET                   , "MET/F");

  m_tree_Mtt->Branch("isSel"              , &m_mtt_isSel                 , "isSel/I");
  m_tree_Mtt->Branch("oneMatchedCombi"    , &m_mtt_OneMatchedCombi       , "oneMatchedCombi/I");
  m_tree_Mtt->Branch("bestSolChi2"        , &m_mtt_BestSolChi2           , "bestSolChi2/F");
  m_tree_Mtt->Branch("isBestSolMatched"   , &m_mtt_IsBestSolMatched      , "isBestSolMatched/I");

  m_tree_Mtt->Branch("numComb"            , &m_mtt_NumComb                , "numComb/I");
  m_tree_Mtt->Branch("solChi2"            , &m_mtt_SolChi2                , "solChi2[numComb]/F");


  m_tree_Mtt->Branch("mLepW_AfterChi2"        , &m_mLepW_AfterChi2       , "mLepW_AfterChi2/F");
  m_tree_Mtt->Branch("mHadW_AfterChi2"        , &m_mHadW_AfterChi2       , "mHadW_AfterChi2/F");

  m_tree_Mtt->Branch("mLepTop_AfterChi2"      , &m_mLepTop_AfterChi2     , "mLepTop_AfterChi2/F");
  m_tree_Mtt->Branch("lepTopPt_AfterChi2"     , &m_lepTopPt_AfterChi2    , "lepTopPt_AfterChi2/F");
  m_tree_Mtt->Branch("lepTopEta_AfterChi2"    , &m_lepTopEta_AfterChi2   , "lepTopEta_AfterChi2/F");
  m_tree_Mtt->Branch("lepTopP4_AfterChi2"     , &m_lepTopP4_AfterChi2);

  m_tree_Mtt->Branch("mHadTop_AfterChi2"      , &m_mHadTop_AfterChi2     , "mHadTop_AfterChi2/F");
  m_tree_Mtt->Branch("hadTopPt_AfterChi2"     , &m_hadTopPt_AfterChi2    , "hadTopPt_AfterChi2/F");
  m_tree_Mtt->Branch("hadTopEta_AfterChi2"    , &m_hadTopEta_AfterChi2   , "hadTopEta_AfterChi2/F");
  m_tree_Mtt->Branch("hadTopP4_AfterChi2"     , &m_hadTopP4_AfterChi2);

  m_tree_Mtt->Branch("pt_tt_AfterChi2"        , &m_pt_tt_AfterChi2       , "pt_tt_AfterChi2/F");
  m_tree_Mtt->Branch("eta_tt_AfterChi2"       , &m_eta_tt_AfterChi2      , "eta_tt_AfterChi2/F");
  m_tree_Mtt->Branch("beta_tt_AfterChi2"      , &m_beta_tt_AfterChi2     , "eta_tt_AfterChi2/F");
  m_tree_Mtt->Branch("mtt_AfterChi2"          , &m_mtt_AfterChi2         , "mtt_AfterChi2/F");

  // Index of selected particles inside respective collection for mtt computation
  m_tree_Mtt->Branch("selectedLeptonIndex"        , &m_selectedLeptonIndex       , "selectedLeptonIndex/I");
  m_tree_Mtt->Branch("selectedLeptonicBIndex"     , &m_selectedLeptonicBIndex    , "selectedLeptonicBIndex/I");
  m_tree_Mtt->Branch("selectedHadronicBIndex"     , &m_selectedHadronicBIndex    , "selectedHadronicBIndex/I");
  m_tree_Mtt->Branch("selectedHadronicFirstJetIndex"  , &m_selectedHadronicFirstJetIndex  , "selectedHadronicFirstJetIndex/I");
  m_tree_Mtt->Branch("selectedHadronicSecondJetIndex" , &m_selectedHadronicSecondJetIndex , "selectedHadronicSecondJetIndex/I");


  m_tree_Mtt->Branch("trigger_passed", &m_trigger_passed, "trigger_passed/O");

  // Weights and errors from differents scale factors
  m_tree_Mtt->Branch("weight", &m_weight, "weight/F");
  m_tree_Mtt->Branch("weight_error_low", &m_weight_error_low, "weight_error_low/F");
  m_tree_Mtt->Branch("weight_error_high", &m_weight_error_high, "weight_error_high/F");

  // Neutrino Pz calculation study
  m_tree_Mtt->Branch("is_neutrino_pz_corrected", &m_is_neutrino_pz_corrected, "is_neutrino_pz_corrected/O");

  // cuts
  m_tree_Mtt->Branch("pass_vertex_cut", &m_pass_vertex_cut, "pass_vertex_cut/I");
  m_tree_Mtt->Branch("pass_met_cut", &m_pass_met_cut, "pass_met_cut/I");
  m_tree_Mtt->Branch("pass_lepton_cut", &m_pass_lepton_cut, "pass_lepton_cut/I");
  m_tree_Mtt->Branch("pass_jet_cut", &m_pass_jet_cut, "pass_jet_cut/I");

  m_MAIN_doSemiMu = cmsswSettings.getParameter<bool>("do_semimu");

  // METSel()
  m_MET_Pt_Min = cmsswSettings.getParameter<edm::ParameterSet>("met").getParameter<double>("pt_min");

  // MuonSel()
  m_MU_Pt_min_loose = cmsswSettings.getParameter<edm::ParameterSet>("muons_loose").getParameter<double>("pt_min");
  m_MU_Eta_max_loose = cmsswSettings.getParameter<edm::ParameterSet>("muons_loose").getParameter<double>("eta_max");
  m_MU_Iso_max_loose = cmsswSettings.getParameter<edm::ParameterSet>("muons_loose").getParameter<double>("isolation_max");

  m_MU_Pt_min = cmsswSettings.getParameter<edm::ParameterSet>("muons_tight").getParameter<double>("pt_min");
  m_MU_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("muons_tight").getParameter<double>("eta_max");
  m_MU_Iso_max = cmsswSettings.getParameter<edm::ParameterSet>("muons_tight").getParameter<double>("isolation_max");

  // ElectronSel()
  m_ELE_Pt_min_loose = cmsswSettings.getParameter<edm::ParameterSet>("electrons_loose").getParameter<double>("pt_min");
  m_ELE_Eta_max_loose = cmsswSettings.getParameter<edm::ParameterSet>("electrons_loose").getParameter<double>("eta_max");
  m_ELE_Iso_max_loose = cmsswSettings.getParameter<edm::ParameterSet>("electrons_loose").getParameter<double>("isolation_max");

  m_ELE_Pt_min = cmsswSettings.getParameter<edm::ParameterSet>("electrons_tight").getParameter<double>("pt_min");
  m_ELE_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("electrons_tight").getParameter<double>("eta_max");
  m_ELE_Iso_max = cmsswSettings.getParameter<edm::ParameterSet>("electrons_tight").getParameter<double>("isolation_max");

  // JetSel()
  m_JET_Pt_min  = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("pt_min");
  m_JET_Eta_max = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("eta_max");
  m_JET_btag_CSVL = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_CSVL");
  m_JET_btag_CSVM = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_CSVM");
  m_JET_btag_CSVT = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_CSVT");
  m_JET_btag_TCHPT = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_TCHPT");
  m_b_tagging_efficiency = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("b_tagging_efficiency");

  m_MAIN_doUseBTag = cmsswSettings.getParameter<edm::ParameterSet>("chi2_sorting").getParameter<bool>("use_btagging");

  std::string fname = "kfparams_semilept.dat";
  m_KinFit = new KinFit(fname, cmsswSettings);

  m_weight = 1.;
}

mtt_analysis::~mtt_analysis()
{
  delete m_KinFit;

  delete m_lepTopP4_AfterChi2;
  delete m_hadTopP4_AfterChi2;
}


//
// First you start with the different physics objects selections
//

// Vertices
int mtt_analysis::VertexSel()
{
  int n_vtx = m_vertex->getSize();

  if (!n_vtx)
    return 2;

  for (int i = 0; i < n_vtx; ++i)
  {
    if (m_vertex->getVtxIsFake(i)) continue;
    if (m_vertex->getVtxNdof(i) < 4) continue;

    return 1;
  }

  return 2;
}


// MET
int mtt_analysis::METSel()
{
  m_mtt_MET = m_jetMet->getMETLorentzVector(0)->Pt();
  if (m_mtt_MET > m_MET_Pt_Min) {
    return 1;
  }

  return 3;
}


int mtt_analysis::MuonSel()
{
  int goodmuidx = -1;

  int n_mu = m_muon->getSize();

  if (!n_mu)
    return 4;

  for (int i = 0; i < n_mu; i++)
  {

    // Selection:
    // - exactly 1 isolated muon
    // - loose muon veto
    // - loose electron veto


    // Selection from https://twiki.cern.ch/twiki/bin/view/CMS/TWikiTopRefEventSel#Muons
    if (! m_muon->getMuisGlobal(i))
      continue;

    TLorentzVector *muP = m_muon->getMuLorentzVector(i);

    if (fabs(muP->Pt()) <= 26)
      continue;

    if (fabs(muP->Eta()) >= 2.1)
      continue;

    if (m_muon->getMunormChi2(i) >= 10.)
      continue;

    if (m_muon->getTrackerLayersWithMeasurements(i) <= 5)
      continue;

    if (m_muon->getGlobalTrackNumberOfValidMuonHits(i) <= 0)
      continue;

    if (m_muon->getNumberOfMatchedStations(i) <= 1)
      continue;

    if (m_muon->getMudB(i) >= 0.2)
      continue;

    if (m_muon->getdZ(i) >= 0.5)
      continue;

    if (m_muon->getMunValPixelHits(i) <= 0)
      continue;

    /* Isolation cut is done in P2PAT. No need to check that */

    m_mtt_MuonPt[m_mtt_NGoodMuons]   = muP->Pt();
    m_mtt_MuonEta[m_mtt_NGoodMuons]  = muP->Eta();
    m_mtt_MuRelIso[m_mtt_NGoodMuons] = m_muon->getDeltaBetaCorrectedRelativeIsolation(i);

    ++m_mtt_NGoodMuons;
    goodmuidx = i;
  }

  // No muons
  if (m_mtt_NGoodMuons == 0)
    return 4;

  /*
  if (m_mtt_NGoodMuons == 1 && m_mtt_NLooseGoodMuons > 1)
    return 6;
  */

  // We want only 1 isolated muon
  if (m_mtt_NGoodMuons > 1)
    return 5;

  TLorentzVector* selected_p4 = m_muon->getMuLorentzVector(goodmuidx);

  // Muon veto
  int size = m_muon_loose->getSize();
  for (int i = 0; i < size; i++) {

    // Exclude the isolated muon
    TLorentzVector* p4 = m_muon_loose->getMuLorentzVector(i);
    if (*selected_p4 == *p4) {
      continue;
    }

    if ((m_muon_loose->getMuisGlobal(i) || m_muon_loose->getMuisTracker(i)) &&
        (p4->Pt() > 10) &&
        (p4->Eta() < 2.5) &&
        (m_muon_loose->getDeltaBetaCorrectedRelativeIsolation(i) < 0.20)) {

      m_mtt_NLooseGoodMuons++;
    }
  }
  m_mtt_NLooseGoodMuons++; // Our isolated muon is loose too

  if (m_mtt_NLooseGoodMuons > 1)
    return 6;

  //electron veto for semimu channel
  int n_ele = m_electron_loose->getSize();

  for (int i = 0; i < n_ele; ++i) {
    TLorentzVector *eP = m_electron_loose->getEleLorentzVector(i);

    if ((eP->Pt() > 20) &&
        (eP->Eta() < 2.5) &&
        (m_electron_loose->passVetoId(i)) &&
        (m_electron_loose->getRhoCorrectedRelativeIsolation(i) < 0.15))
      return 7;
  }

  m_refLept = m_muon->getMuLorentzVector(goodmuidx);
  m_selectedLeptonIndex = goodmuidx;

  if (m_isMC) {
    // Get scale factor
    ScaleFactor sf = m_muon->getScaleFactor(goodmuidx);
    m_weight *= sf.getValue();
    m_weight_error_low += sf.getErrorLow() * sf.getErrorLow();
    m_weight_error_high += sf.getErrorHigh() * sf.getErrorHigh();
  }

  return 1;
}

int mtt_analysis::ElectronSel()
{
  int goodelidx = -1;

  int n_ele = m_electron->getSize();

  if (!n_ele) return 4;

  for (int i = 0; i < n_ele; i++)
  {
    TLorentzVector *eP = m_electron->getEleLorentzVector(i);

    if (fabs(eP->Pt()) <= 30)
      continue;

    if (fabs(eP->Eta()) >= 2.5)
      continue;

    if (fabs(m_electron->getSuperClusterEta(i)) >= 1.4442 && fabs(m_electron->getSuperClusterEta(i)) < 1.5660)
      continue;

    if (! m_electron->passTightId(i))
      continue;

    m_mtt_ElectronPt[m_mtt_NGoodElectrons] = eP->Pt();
    m_mtt_ElectronEta[m_mtt_NGoodElectrons] = eP->Eta();
    m_mtt_ElRelIso[m_mtt_NGoodElectrons]   = m_electron->getRhoCorrectedRelativeIsolation(i);

    m_mtt_NGoodElectrons++;

    goodelidx = i;
  }

  // No electrons? bye bye
  if (m_mtt_NGoodElectrons == 0)
    return 4;

  // Only one good electron per event
  if (m_mtt_NGoodElectrons > 1)
    return 5;

  TLorentzVector* selected_p4 = m_electron->getEleLorentzVector(goodelidx);

  // Electron veto
  int size = m_electron_loose->getSize();
  for (int i = 0; i < size; i++) {
    TLorentzVector *eP = m_electron_loose->getEleLorentzVector(i);

    if (*selected_p4 == *eP)
      continue;

    if ((eP->Pt() > 20) &&
        (eP->Eta() < 2.5) &&
        (m_electron_loose->passVetoId(i)) &&
        (m_electron_loose->getRhoCorrectedRelativeIsolation(i) < 0.15))
      return 7;
  }

  // Muon veto
  size = m_muon_loose->getSize();
  for (int i = 0; i < size; i++) {
    TLorentzVector* p4 = m_muon_loose->getMuLorentzVector(i);

    if ((m_muon_loose->getMuisGlobal(i) || m_muon_loose->getMuisTracker(i)) &&
        (p4->Pt() > 10) &&
        (p4->Eta() < 2.5) &&
        (m_muon_loose->getDeltaBetaCorrectedRelativeIsolation(i) < 0.20)) {
      return 7;
    }
  }

  m_refLept = m_electron->getEleLorentzVector(goodelidx);
  m_selectedLeptonIndex = goodelidx;

  if (m_isMC) {
    // Get scale factor
    ScaleFactor sf = m_electron->getScaleFactor(goodelidx);
    m_weight *= sf.getValue();
    m_weight_error_low += sf.getErrorLow() * sf.getErrorLow();
    m_weight_error_high += sf.getErrorHigh() * sf.getErrorHigh();
  }

  return 1;
}


int mtt_analysis::JetSel()
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
    TLorentzVector *jetP = m_jetMet->getJetLorentzVector(i);

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

#define CHECK_RES_AND_RETURN(res, var) \
  if (res != 1) { \
    var = 0; \
    m_mtt_isSel = res; \
    fillTree(); \
    return; \
  } else { \
    var = 1; \
  }

#define   MTT_TRIGGER_NOT_FOUND   1000

void mtt_analysis::analyze(const edm::Event& event, const edm::EventSetup& iSetup, PatExtractor& extractor) {
  analyze(iSetup, extractor);
}

void mtt_analysis::analyze(const edm::EventSetup& iSetup, PatExtractor& extractor)
{
  reset();

  m_refLept  = nullptr;

  m_vertex   = std::static_pointer_cast<VertexExtractor>(extractor.getExtractor("vertex"));
  
  m_muon     = std::static_pointer_cast<MuonExtractor>(extractor.getExtractor("muons"));
  m_muon_loose = std::static_pointer_cast<MuonExtractor>(extractor.getExtractor("muons_loose"));

  m_electron = std::static_pointer_cast<ElectronExtractor>(extractor.getExtractor("electrons"));
  m_electron_loose = std::static_pointer_cast<ElectronExtractor>(extractor.getExtractor("electrons_loose"));

  m_jetMet   = std::static_pointer_cast<JetMETExtractor>(extractor.getExtractor("JetMET"));

  m_event    = std::static_pointer_cast<EventExtractor>(extractor.getExtractor("event"));

  std::shared_ptr<HLTExtractor> HLT = std::static_pointer_cast<HLTExtractor>(extractor.getExtractor("HLT"));
  m_trigger_passed = HLT->isTriggerFired();

  if (m_isMC)
  {
    m_MC = std::static_pointer_cast<MCExtractor>(extractor.getExtractor("MC"));
    MCidentification();
  }

  m_nPU = m_event->nPU();

  int res = VertexSel();
  CHECK_RES_AND_RETURN(res, m_pass_vertex_cut);

  res = METSel();
  CHECK_RES_AND_RETURN(res, m_pass_met_cut);

  if (m_MAIN_doSemiMu) {
    res  = MuonSel();
  } else {
    res = ElectronSel();
  }
  CHECK_RES_AND_RETURN(res, m_pass_lepton_cut);

  res = JetSel();
  CHECK_RES_AND_RETURN(res, m_pass_jet_cut);

  m_mtt_isSel = 1;

  // We selected a candidate (and m_refLept)
  loopOverCombinations();

  fillTree();
}



void mtt_analysis::loopOverCombinations()
{
  //jets indices
  int c_j1 = -1;
  int c_j2 = -1;
  int c_j3 = -1;
  int c_j4 = -1;

  //
  int bestj1 = -1;
  int bestj2 = -1;
  int bestj3 = -1;
  int bestj4 = -1;

  int n_btaggedjets = 0;

  //chi2 variables
  double minfitchi2 = std::numeric_limits<double>::infinity();
  double fitchi2    = std::numeric_limits<double>::infinity();

  m_mtt_NumComb = 0;

  int n_jets = m_selJetsIds.size();

  //If requested get the number of b-tagged jets in the selected jets sample
  if (m_MAIN_doUseBTag)
    n_btaggedjets = m_mtt_NBtaggedJets_CSVM;

  //std::cout << "Nulber of b tagged jets: " << btaggedjets.size() << std::endl;

  /**
   * Indices:
   * Jet 1: First b-jet (hadronic)
   * Jet 2: Second b-jet (leptonic)
   * Jet 3: First light jet
   * Jet 4: Second light jet
   */

  //calculate all the possible jets combinations for chi2 calculation.
  //the idea is: b-tag fake rate is low, so b-tagged jets are only associated to bjets (indices bjet1idx and bjet2idx below)
  //on the other hand the b-tag efficiency is low, so non-b-tagged jets can be associated to bjets
  //bottom line: only non b-tagged jets can be associated to the light jets (indices jet3idx and jet4idx below)
  //in the chi2 calculation for this reason in the selected jets we need at least 2 non b-tagged jets

  int numberoflightjets = n_jets - n_btaggedjets;

  if (numberoflightjets < 2)
    return; // if we dont have at least 2 non b-tagged jets, chi2 is -1

  for (int bj1 = 0; bj1 < n_jets; ++bj1)
  {
    for (int bj2 = 0; bj2 < n_jets; ++bj2)
    {
      if (bj2 == bj1)
        continue; //dont pick the one you already used

      for (int j3 = 0; j3 < n_jets; ++j3)
      {
        // dont pick the two jets you used, or btagged jets
        if (j3 == bj1 || j3 == bj2 || (m_MAIN_doUseBTag && isBJet(m_selJetsIds[j3])))
          continue;

        for (int j4 = j3 + 1; j4 < n_jets; ++j4)
        {
          //dont pick the two jets you used, or btagged jets
          if (j4 == bj1 || j4 == bj2 || (m_MAIN_doUseBTag && isBJet(m_selJetsIds[j4])))
            continue;

          // Get the real jet indices

          c_j1 = m_selJetsIds[bj1];
          c_j2 = m_selJetsIds[bj2];
          c_j3 = m_selJetsIds[j3];
          c_j4 = m_selJetsIds[j4];

          /// Try to find a matching solution (no lept for the moment)

          if (m_isMC)
            m_mtt_OneMatchedCombi = match_MC(c_j1, c_j2, c_j3, c_j4, 0);

          // This call corrects MET pz
          bool res = m_KinFit->ReadObjects(*m_jetMet->getJetLorentzVector(c_j3),
              *m_jetMet->getJetLorentzVector(c_j4),
              *m_jetMet->getJetLorentzVector(c_j1),
              *m_refLept,
              *m_jetMet->getMETLorentzVector(0),
              *m_jetMet->getJetLorentzVector(c_j2),
              m_MAIN_doSemiMu
              );

          if (!res)
            return; // We will never get anything with this event

          /*
          if (m_MAIN_doKF) //use the kinfit to choose the best jet pairing
          {
            (m_KinFit->Fit()) // Do the kinfit converged
              ? fitchi2 = m_KinFit->GetKFChi2()
              : fitchi2 = std::numeric_limits<double>::infinity();
          }
          else*/  //else use the chi2 to chose the best combination and apply the kinfit only to this latter
          {
            fitchi2 = m_KinFit->GlobalSimpleChi2(AllJetsPt);
            //std::cout << "Chi2: " << fitchi2 << std::endl;
          }

          if (fitchi2 < minfitchi2)
          {
            minfitchi2 = fitchi2;
            bestj1     = c_j1;
            bestj2     = c_j2;
            bestj3     = c_j3;
            bestj4     = c_j4;

            m_mtt_SolChi2[m_mtt_NumComb] = fitchi2;
            ++m_mtt_NumComb;
          }

        } // j4 close
      } // j3 close
    } // j2 close
  } // we out of the loop over jet pairings

  m_selectedHadronicBIndex = bestj1;
  m_selectedLeptonicBIndex = bestj2;
  m_selectedHadronicFirstJetIndex = bestj3;
  m_selectedHadronicSecondJetIndex = bestj4;

  // Do second KinFit
  // We looped over the combinations and do the kinfit (if not already done)
  if (m_mtt_NumComb > 0)
  {
    
    // Put selected object inside KinFit. This will correct MET and everything we need
    m_KinFit->ReadObjects(*m_jetMet->getJetLorentzVector(bestj3),
        *m_jetMet->getJetLorentzVector(bestj4),
        *m_jetMet->getJetLorentzVector(bestj1),
        *m_refLept,
        *m_jetMet->getMETLorentzVector(0),
        *m_jetMet->getJetLorentzVector(bestj2),
        m_MAIN_doSemiMu,
        &m_is_neutrino_pz_corrected
        );

    const TLorentzVector& measuredLepton = m_KinFit->GetMeasuredLepton();
    const TLorentzVector& measuredNeutrino = m_KinFit->GetMeasuredNeutrino();
    const TLorentzVector& measuredLeptonicB = m_KinFit->GetMeasuredLeptonicBJet();
    const TLorentzVector& measuredHadronicB = m_KinFit->GetMeasuredHadronicBJet();
    const TLorentzVector& measuredHadronicFirstJet = m_KinFit->GetMeasuredFirstLightJet();
    const TLorentzVector& measuredHadronicSecondJet = m_KinFit->GetMeasuredSecondLightJet();

    /**
     * Compute Mtt before doing KinFit
     */
    m_mLepW_AfterChi2   = (measuredNeutrino + measuredLepton).M();
    m_mHadW_AfterChi2   = (measuredHadronicFirstJet + measuredHadronicSecondJet).M();

    *m_lepTopP4_AfterChi2 = (measuredLepton + measuredNeutrino + measuredLeptonicB);
    m_mLepTop_AfterChi2 = m_lepTopP4_AfterChi2->M();
    m_lepTopPt_AfterChi2 = m_lepTopP4_AfterChi2->Pt();
    m_lepTopEta_AfterChi2 = m_lepTopP4_AfterChi2->Eta();

    *m_hadTopP4_AfterChi2 = (measuredHadronicFirstJet + measuredHadronicSecondJet + measuredHadronicB);
    m_mHadTop_AfterChi2 = m_hadTopP4_AfterChi2->M();
    m_hadTopPt_AfterChi2 = m_hadTopP4_AfterChi2->Pt();
    m_hadTopEta_AfterChi2 = m_hadTopP4_AfterChi2->Eta();

    TLorentzVector res = (*m_lepTopP4_AfterChi2 + *m_hadTopP4_AfterChi2);
    m_mtt_AfterChi2     = res.M();
    m_pt_tt_AfterChi2   = res.Pt();
    m_eta_tt_AfterChi2   = res.Eta();
    m_beta_tt_AfterChi2  = fabs(res.Pz() / res.E());

    /*
       std::cout << "Mt lepton: " << m_mLepTop_AfterChi2 << std::endl;
       std::cout << "Mt hadronic: " << m_mHadTop_AfterChi2 << std::endl;
       std::cout << "Mtt: " << m_mtt_AfterChi2 << std::endl;
       */
    /*

    m_KinFit->ReadObjects(*m_jetMet->getJetLorentzVector(bestj3),
        *m_jetMet->getJetLorentzVector(bestj4),
        *m_jetMet->getJetLorentzVector(bestj1),
        *m_refLept,
        *m_jetMet->getMETLorentzVector(0),
        *m_jetMet->getJetLorentzVector(bestj2),
        m_MAIN_doSemiMu
        );

    (m_KinFit->Fit()) // Do the kinfit converged
      ? fitchi2 = m_KinFit->GetKFChi2()
      : fitchi2 = std::numeric_limits<double>::infinity();

    */
  }

  if (m_isMC)
    m_mtt_IsBestSolMatched = match_MC(bestj1, bestj2, bestj3, bestj4, 0);

  //m_mtt_KFChi2 = fitchi2;

  //if (!m_MAIN_doKF)
  m_mtt_BestSolChi2 = minfitchi2;

  /*
  if (fitchi2 <= 1.E6)
  {
    m_mLepTop_AfterChi2andKF = (*m_KinFit->GetFittedLeptonicBJet() + *m_KinFit->GetFittedLepton() + *m_KinFit->GetFittedNeutrino()).M();
    m_mHadTop_AfterChi2andKF = (*m_KinFit->GetFittedHadronicBJet() + *m_KinFit->GetFittedFirstLightJet() + *m_KinFit->GetFittedSecondLightJet()).M();
    m_mtt_AfterChi2andKF     = (*m_KinFit->GetFittedLeptonicBJet() + *m_KinFit->GetFittedLepton() + *m_KinFit->GetFittedNeutrino() + *m_KinFit->GetFittedHadronicBJet() + *m_KinFit->GetFittedFirstLightJet() + *m_KinFit->GetFittedSecondLightJet()).M();
  }
  else
  {
    m_mLepTop_AfterChi2andKF = std::numeric_limits<double>::infinity();
    m_mHadTop_AfterChi2andKF = std::numeric_limits<double>::infinity();
    m_mtt_AfterChi2andKF     = std::numeric_limits<double>::infinity();
  }
  */
}


// MC stuff

#define ID_B (5)
#define ID_T (6)

#define ID_E (11)
#define ID_NEUTRINO_E (12)
#define ID_MU (13)
#define ID_NEUTRINO_MU (14)
#define ID_TAU (15)
#define ID_NEUTRINO_TAU (16)

#define ID_W (24)

int mtt_analysis::patIndexToExtractorIndex(int patIndex) const {

  for (int i = 0; i < m_MC->getSize() ; i++) {
    if (m_MC->getPatIndex(i) == patIndex)
      return i;
  }

  return -1;
}

void mtt_analysis::MCidentification()
{
  nEle    = 0;
  nMu     = 0;
  nTau    = 0;
  nNuEle  = 0;
  nNuMu   = 0;
  nNuTau  = 0;
  nQuarkb = 0;
  nTop    = 0;
  Top.clear();
  std::vector<int> topIndexes;

  int n_MC = m_MC->getSize();

  if (!n_MC)
    return;

  for (int i = 0; i < n_MC ; ++i)
  {

    if (abs(m_MC->getType(i)) == ID_T) {

      if (std::find(topIndexes.begin(), topIndexes.end(), i) != topIndexes.end())
        continue;

      topIndexes.push_back(i);
      TLorentzVector TL_Top(m_MC->getPx(i),
          m_MC->getPy(i),
          m_MC->getPz(i),
          m_MC->getE(i));

      Top.push_back(TL_Top);
      nTop++;

      continue;
    }


    int motherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(i));
    int grandMotherIndex = -1;
    if (motherIndex != -1)
      grandMotherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(motherIndex));

    if (motherIndex == -1)
      continue;

    if (abs(m_MC->getType(motherIndex)) == ID_T || (grandMotherIndex != -1 && abs(m_MC->getType(grandMotherIndex) == ID_T)))
    {
      // Skip W
      if (abs(m_MC->getType(motherIndex)) == ID_W)
        continue;

      /// Count the number of leptons and neutrinos from Top->W
      if (fabs(m_MC->getType(i)) == 11) ++nEle;   //Electron from WTop
      if (fabs(m_MC->getType(i)) == 13) ++nMu;    //Muon	   from WTop
      if (fabs(m_MC->getType(i)) == 15) ++nTau;   //Tau	   from WTop
      if (fabs(m_MC->getType(i)) == 12) ++nNuEle; //NuEle    from WTop
      if (fabs(m_MC->getType(i)) == 14) ++nNuMu;  //NuMu	   from WTop
      if (fabs(m_MC->getType(i)) == 16) ++nNuTau; //NuTau    from WTop
    }

    /// Count the number of b quark from Top
    if (abs(m_MC->getType(i)) == ID_B && abs(m_MC->getType(motherIndex)) == ID_T)
    {
      nQuarkb++; //Quark b from Top
    }
  }

  if (nEle == 1 && nNuEle == 1 && nMu == 0 && nNuMu == 0 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 1;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 1 && nNuMu == 1 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 2;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 0 && nNuMu == 0 && nTau == 1 && nNuTau == 1 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 3;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 0 && nNuMu == 0 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 4;
  }

  if (nEle == 2 && nNuEle == 2 && nMu == 0 && nNuMu == 0 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 5;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 2 && nNuMu == 2 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 6;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 0 && nNuMu == 0 && nTau == 2 && nNuTau == 2 && nQuarkb > 1 && nTop == 2)
  {
    m_MC_channel = 7;
  }

  if (nEle == 1 && nNuEle == 1 && nMu == 1 && nNuMu == 1 && nTau == 0 && nNuTau == 0 && nQuarkb == 2 && nTop == 2)
  {
    m_MC_channel = 8 ;
  }

  if (nEle == 1 && nNuEle == 1 && nMu == 0 && nNuMu == 0 && nTau == 1 && nNuTau == 1 && nQuarkb == 2 && nTop == 2)
  {
    m_MC_channel = 9 ;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 1 && nNuMu == 1 && nTau == 1 && nNuTau == 1 && nQuarkb == 2 && nTop == 2)
  {
    m_MC_channel = 10;
  }

  if (nTop == 2) {

    TLorentzVector mc_resonance = Top[0] + Top[1];

    m_MC_mtt = mc_resonance.M();
    m_MC_pt_tt = mc_resonance.Pt();
    m_MC_eta_tt = mc_resonance.Eta();
    m_MC_beta_tt = fabs(mc_resonance.Pz() / mc_resonance.E());
  }

  if (m_MC_channel != 1 && m_MC_channel != 2) {
    return;
  }

  // Extract index of semi-leptonic event, and store them in tree. Useful if you want to know how many jets you have selected right
  if (false) {
    std::cout << "New event" << std::endl;
    for (int i = 0; i < n_MC; i++) {
      std::cout << "\t[" << i << "] Type: " << m_MC->getType(i) << std::endl;
    }
  }

  bool keepEvent = true;
  for (int i = 0; i < n_MC; i++) {
    
    int motherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(i));
    int grandMotherIndex = -1;
    if (motherIndex != -1)
      grandMotherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(motherIndex));

    if (motherIndex == -1 || grandMotherIndex == -1)
      continue;

    // Look only event coming (directly / indirectly) from a top
    if (abs(m_MC->getType(motherIndex)) == ID_T || abs(m_MC->getType(grandMotherIndex)) == ID_T)  {

      int type = abs(m_MC->getType(i));
      // W? Continue
      if (type == ID_W)
        continue;

      // Only semi-mu or semi-e events are interesting, so throw away event with a tau
      if (type == ID_TAU) {
        keepEvent = false;
        break;
      }

      switch (type) {
        case ID_E:
          if (m_leptonIndex != -1) {
            keepEvent = false;
            break;
          }
          m_leptonIndex = i;
          break;

        case ID_MU:
          if (m_leptonIndex != -1) {
            keepEvent = false;
            break;
          }
          m_leptonIndex = i;
          break;

        case ID_NEUTRINO_E:
        case ID_NEUTRINO_MU:
        case ID_NEUTRINO_TAU:
          if (m_neutrinoIndex != -1) {
            keepEvent = false;
            break;
          }
          m_neutrinoIndex = i;

          // In some case, Madgraph does not output a W in the LHE
          if (abs(m_MC->getType(grandMotherIndex)) == ID_T)
            m_leptonicTopIndex = grandMotherIndex;
          else
            m_leptonicTopIndex = motherIndex;

          break;

        case ID_B:
          if (m_leptonicBIndex == -1) {
            m_leptonicBIndex = i;
          } else {
            if (m_hadronicBIndex != -1) {
              keepEvent = false;
              break;
            }
            m_hadronicBIndex = i;
          }
          break;

        default: // Other jets
          if (m_firstJetIndex == -1) {
            m_firstJetIndex = i;
          } else {
            if (m_secondJetIndex != -1) {
              keepEvent = false;
              break;
            }
            m_secondJetIndex = i;
          }
          break;
      }

      if (! keepEvent)
        break;
    }
  }

  if (m_leptonIndex == -1 || m_neutrinoIndex == -1 || m_leptonicBIndex == -1 || m_hadronicBIndex == -1)
    keepEvent = false;

  if (! keepEvent) {
    m_leptonIndex = m_leptonicBIndex = m_hadronicBIndex = m_neutrinoIndex = m_firstJetIndex = m_secondJetIndex = m_leptonicTopIndex = -1;
    return;
  }

  // Reorder B jet indexes
  if (patIndexToExtractorIndex(m_MC->getMom1Index(m_leptonicBIndex)) != m_leptonicTopIndex) {
    // Wrong combinaison, swap
    std::swap(m_leptonicBIndex, m_hadronicBIndex);
  }

  if (false) {
    std::cout << "Lepton index: " << m_leptonIndex << std::endl;
    std::cout << "Neutrino index: " << m_neutrinoIndex << std::endl;
    std::cout << "Leptonic B index: " << m_leptonicBIndex << std::endl;
    std::cout << "Hadronic B index: " << m_hadronicBIndex << std::endl;
    std::cout << "First jet index: " << m_firstJetIndex << std::endl;
    std::cout << "Second jet index: " << m_secondJetIndex << std::endl;
  }

  // Compute masses
  TLorentzVector mc_hadr_W = *m_MC->getP4(m_firstJetIndex) + *m_MC->getP4(m_secondJetIndex);
  m_MC_hadronicWMass = mc_hadr_W.M();

  TLorentzVector mc_hadr_top = mc_hadr_W + *m_MC->getP4(m_hadronicBIndex);
  m_MC_hadronicTopMass = mc_hadr_top.M();

  TLorentzVector mc_lept_W = *m_MC->getP4(m_neutrinoIndex) + *m_MC->getP4(m_leptonIndex);
  m_MC_leptonicWMass = mc_lept_W.M();

  TLorentzVector mc_lept_top = mc_lept_W + *m_MC->getP4(m_leptonicBIndex);
  m_MC_leptonicTopMass = mc_lept_top.M();
}

int mtt_analysis::match_MC(int idxJetbH, int idxJetbL, int idxJet1,	int idxJet2,
    int idxLepton)
{
  if (
      /// Ask if Jet b hadronique  come from a b and top
      fabs(m_MC->getType(m_jetMet->getJetMCIndex(idxJetbH))) == 5 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_jetMet->getJetMCIndex(idxJetbH)))) == 6 &&

      /// Ask if Jet b leptonique  come from a b and top
      fabs(m_MC->getType(m_jetMet->getJetMCIndex(idxJetbL))) == 5 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_jetMet->getJetMCIndex(idxJetbL)))) == 6 &&

      /// Ask if jet 1,2 come from light quark and W and top
      fabs(m_MC->getType(m_jetMet->getJetMCIndex(idxJet1))) < 5 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_jetMet->getJetMCIndex(idxJet1)))) == 24 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_MC->getMom1Index(m_jetMet->getJetMCIndex(idxJet1))))) == 6 &&
      fabs(m_MC->getType(m_jetMet->getJetMCIndex(idxJet2))) < 5 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_jetMet->getJetMCIndex(idxJet2)))) == 24 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_MC->getMom1Index(m_jetMet->getJetMCIndex(idxJet2))))) == 6
     )
  {
    return 1;
  }
  else
  {
    return 0;
  }
}


void mtt_analysis::fillTree()
{
  m_weight_error_low = sqrt(m_weight_error_low);
  m_weight_error_high = sqrt(m_weight_error_high);
  m_tree_Mtt->Fill();
}

// Here we just reset the ROOTtree parameters

void mtt_analysis::reset()
{
  m_pass_vertex_cut = -1;
  m_pass_met_cut = -1;
  m_pass_lepton_cut = -1;
  m_pass_jet_cut = -1;

  m_mtt_isSel = 0;
  m_mtt_IsBestSolMatched = -1;
  m_mtt_OneMatchedCombi = 0;
  m_mtt_BestSolChi2 = -1.;
  //m_mtt_KFChi2 = -1.;
  //m_mtt_AfterChi2andKF = -1.;
  m_mtt_NumComb = 0;
  //m_mLepTop_AfterChi2andKF = -1.;
  //m_mHadTop_AfterChi2andKF = -1.;

  m_mtt_AfterChi2     = -1.;
  m_pt_tt_AfterChi2   = -1.;
  m_eta_tt_AfterChi2  = -1.;
  m_beta_tt_AfterChi2 = -1.;
  m_mLepTop_AfterChi2 = -1.;
  m_mHadTop_AfterChi2 = -1.;
  m_mHadW_AfterChi2   = -1.;

  m_lepTopPt_AfterChi2  = -1;
  m_lepTopEta_AfterChi2 = -1;
  m_hadTopPt_AfterChi2  = -1;
  m_hadTopEta_AfterChi2 = -1;

  m_lepTopP4_AfterChi2->SetPxPyPzE(0., 0., 0., 0.);
  m_hadTopP4_AfterChi2->SetPxPyPzE(0., 0., 0., 0.);

  m_selectedLeptonIndex               = -1;
  m_selectedLeptonicBIndex            = -1;
  m_selectedHadronicBIndex            = -1;
  m_selectedHadronicFirstJetIndex     = -1;
  m_selectedHadronicSecondJetIndex    = -1;

  m_mtt_NGoodMuons = 0;
  m_mtt_NLooseGoodMuons = 0;
  m_mtt_NGoodElectrons = 0;

  for (int tmp = 0; tmp < 20; ++tmp)
  {
    m_mtt_MuonPt[tmp] = 0.;
    m_mtt_2DDrMin[tmp] = -1.;
    m_mtt_2DpTrel[tmp] = -1.;
    m_mtt_MuRelIso[tmp] = -1.;

    m_mtt_ElectronPt[tmp] = 0.;
    m_mtt_HyperTight1MC[tmp] = -1;
    m_mtt_ElRelIso[tmp] = -1.;
  }


  m_mtt_NJets = 0;
  m_mtt_NBtaggedJets_CSVL = 0;
  m_mtt_NBtaggedJets_CSVM = 0;
  m_mtt_NBtaggedJets_CSVT = 0;
  //m_mtt_NBtaggedJets_TCHPL = 0;
  //m_mtt_NBtaggedJets_TCHPM = 0;
  m_mtt_NBtaggedJets_TCHPT = 0;
  //m_mtt_NBtaggedJets_SSVHEM = 0;
  //m_mtt_NBtaggedJets_SSVHPT = 0;

  m_mtt_1stjetpt = 0.;
  m_mtt_2ndjetpt = 0.;
  m_mtt_3rdjetpt = 0.;
  m_mtt_4thjetpt = 0.;

  for (int tmp = 0; tmp < 100; ++tmp)
  {
    m_mtt_JetEta[tmp]     = 1000.;
    m_mtt_JetPt[tmp]      = 0.;
  }

  m_mtt_MET = 0.;
  m_nPU        = 0.;
  m_MC_channel = 0;
  m_MC_mtt     = -1.;

  m_leptonIndex = -1;
  m_neutrinoIndex = -1;

  m_leptonicBIndex = -1;
  m_hadronicBIndex = -1;
  m_leptonicTopIndex = -1;

  m_firstJetIndex = -1;
  m_secondJetIndex = -1;
  
  m_trigger_passed = false;

  m_MC_pt_tt         = -1;
  m_MC_eta_tt        = -1;
  m_MC_beta_tt       = -1;
  m_MC_hadronicWMass = -1;
  m_MC_leptonicWMass = -1;
  m_MC_hadronicTopMass = -1;
  m_MC_leptonicTopMass = -1;

  m_weight = 1.;
  m_weight_error_low = 0.;
  m_weight_error_high = 0.;

  m_is_neutrino_pz_corrected = false;
}

bool mtt_analysis::isBJet(unsigned int index) {
  // Use recommanded WP from https://indico.cern.ch/getFile.py/access?contribId=4&resId=2&materialId=slides&confId=195042
  return m_jetMet->getJetBTagProb_CSV(index) > m_JET_btag_CSVM;
}

}

DEFINE_EDM_PLUGIN(PatExtractorPluginFactory, patextractor::mtt_analysis, "mtt_analysis");
