#include <iostream>
#include <vector>
#include <limits>

#include "TH2.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "Extractors/PatExtractor/interface/AnalysisSettings.h"
#include "Extractors/PatExtractor/interface/MCExtractor.h"
#include "Extractors/PatExtractor/interface/HLTExtractor.h"
#include "Extractors/PatExtractor/interface/MuonExtractor.h"
#include "Extractors/PatExtractor/interface/ElectronExtractor.h"
#include "Extractors/PatExtractor/interface/JetExtractor.h"
#include "Extractors/PatExtractor/interface/METExtractor.h"
#include "Extractors/PatExtractor/interface/VertexExtractor.h"
#include "Extractors/PatExtractor/interface/KinFit.h"
#include "Extractors/PatExtractor/interface/EventExtractor.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include <FWCore/Framework/interface/EventSetup.h>

#include "../interface/mtt_analysis_new.h"

using namespace std;

mtt_analysis_new::mtt_analysis_new(AnalysisSettings *settings)
{
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
  m_tree_Mtt->Branch("MC_leptonicWIndex"  , &m_leptonicWIndex         , "MC_leptonicWIndex/I");
  m_tree_Mtt->Branch("MC_leptonicBIndex"  , &m_leptonicBIndex         , "MC_leptonicBIndex/I");

  m_tree_Mtt->Branch("MC_hadronicBIndex"  , &m_hadronicBIndex         , "MC_hadronicBIndex/I");
  m_tree_Mtt->Branch("MC_hadronicFirstJetIndex" , &m_firstJetIndex    , "MC_hadronicFirstJetIndex/I");
  m_tree_Mtt->Branch("MC_hadronicSecondJetIndex", &m_secondJetIndex   , "MC_hadronicSecondJetIndex/I");

  m_tree_Mtt->Branch("nGoodMuons"         , &m_mtt_NGoodMuons         , "nGoodMuons/I");
  m_tree_Mtt->Branch("nLooseGoodMuons"    , &m_mtt_NLooseGoodMuons    , "nLooseGoodMuons/I");
  m_tree_Mtt->Branch("muonPt"             , &m_mtt_MuonPt             , "muonPt[nGoodMuons]/F");
  m_tree_Mtt->Branch("2DDrMin"            , &m_mtt_2DDrMin            , "2DDrMin[nGoodMuons]/F");
  m_tree_Mtt->Branch("2DpTrel"            , &m_mtt_2DpTrel            , "2DpTrel[nGoodMuons]/F");
  m_tree_Mtt->Branch("muRelIso"           , &m_mtt_MuRelIso           , "muRelIso[nGoodMuons]/F");

  m_tree_Mtt->Branch("nGoodElectrons"     , &m_mtt_NGoodElectrons     , "nGoodElectrons/I");
  m_tree_Mtt->Branch("electronPt"         , &m_mtt_ElectronPt         , "electronPt[nGoodElectrons]/F");
  m_tree_Mtt->Branch("elRelIso"           , &m_mtt_ElRelIso           , "elRelIso[nGoodElectrons]/F");
  m_tree_Mtt->Branch("hyperTight1MC"      , &m_mtt_HyperTight1MC      , "hyperTight1MC[nGoodElectrons]/I");


  m_tree_Mtt->Branch("1stjetpt"           , &m_mtt_1stjetpt           , "1stjetpt/F");
  m_tree_Mtt->Branch("2ndjetpt"           , &m_mtt_2ndjetpt           , "2ndjetpt/F");
  m_tree_Mtt->Branch("3rdjetpt"           , &m_mtt_3rdjetpt           , "3rdjetpt/F");
  m_tree_Mtt->Branch("4thjetpt"           , &m_mtt_4thjetpt           , "4thjetpt/F");
  m_tree_Mtt->Branch("nJets"              , &m_mtt_NJets              , "nJets/I");
  m_tree_Mtt->Branch("jetEta"             , &m_mtt_JetEta             , "jetEta[nJets]/F");
  m_tree_Mtt->Branch("jetPt"              , &m_mtt_JetPt              , "jetPt[nJets]/F");

  //m_tree_Mtt->Branch("nBtaggedJets_TCHEL" , &m_mtt_NBtaggedJets_TCHEL    , "nBtaggedJets_TCHEL/I");
  //m_tree_Mtt->Branch("nBtaggedJets_TCHEM" , &m_mtt_NBtaggedJets_TCHEM    , "nBtaggedJets_TCHEM/I");
  //m_tree_Mtt->Branch("nBtaggedJets_TCHET" , &m_mtt_NBtaggedJets_TCHET    , "nBtaggedJets_TCHET/I");
  //m_tree_Mtt->Branch("nBtaggedJets_TCHPL" , &m_mtt_NBtaggedJets_TCHPL    , "nBtaggedJets_TCHPL/I");
  //m_tree_Mtt->Branch("nBtaggedJets_TCHPM" , &m_mtt_NBtaggedJets_TCHPM    , "nBtaggedJets_TCHPM/I");
  m_tree_Mtt->Branch("nBtaggedJets_TCHPT" , &m_mtt_NBtaggedJets_TCHPT    , "nBtaggedJets_TCHPT/I");
  //m_tree_Mtt->Branch("nBtaggedJets_SSVHEM", &m_mtt_NBtaggedJets_SSVHEM   , "nBtaggedJets_SSVHEM/I");
  //m_tree_Mtt->Branch("nBtaggedJets_SSVHPT", &m_mtt_NBtaggedJets_SSVHPT   , "nBtaggedJets_SSVHPT/I");
  m_tree_Mtt->Branch("nBtaggedJets_CSVL" ,  &m_mtt_NBtaggedJets_CSVL    , "nBtaggedJets_CSVL/I");
  m_tree_Mtt->Branch("nBtaggedJets_CSVM" ,  &m_mtt_NBtaggedJets_CSVM    , "nBtaggedJets_CSVM/I");
  m_tree_Mtt->Branch("nBtaggedJets_CSVT" ,  &m_mtt_NBtaggedJets_CSVT    , "nBtaggedJets_CSVT/I");

  m_tree_Mtt->Branch("MET"                , &m_mtt_MET                   , "MET/F");

  m_tree_Mtt->Branch("isSel"              , &m_mtt_isSel                 , "isSel/I");
  m_tree_Mtt->Branch("oneMatchedCombi"    , &m_mtt_OneMatchedCombi       , "oneMatchedCombi/I");
  m_tree_Mtt->Branch("bestSolChi2"        , &m_mtt_BestSolChi2           , "bestSolChi2/F");
  m_tree_Mtt->Branch("isBestSolMatched"   , &m_mtt_IsBestSolMatched      , "isBestSolMatched/I");
  m_tree_Mtt->Branch("KFChi2"             , &m_mtt_KFChi2                , "KFChi2/F");

  m_tree_Mtt->Branch("numComb"            , &m_mtt_NumComb                , "numComb/I");
  m_tree_Mtt->Branch("solChi2"            , &m_mtt_SolChi2                , "solChi2[numComb]/F");


  m_tree_Mtt->Branch("mHadW_AfterChi2"        , &m_mHadW_AfterChi2       , "mHadW_AfterChi2/F");
  m_tree_Mtt->Branch("mLepTop_AfterChi2"      , &m_mLepTop_AfterChi2     , "mLepTop_AfterChi2/F");
  m_tree_Mtt->Branch("mHadTop_AfterChi2"      , &m_mHadTop_AfterChi2     , "mHadTop_AfterChi2/F");
  m_tree_Mtt->Branch("mtt_AfterChi2"          , &m_mtt_AfterChi2         , "mtt_AfterChi2/F");

  m_tree_Mtt->Branch("mLepTop_AfterChi2andKF" , &m_mLepTop_AfterChi2andKF, "mLepTop_AfterChi2andKF/F");
  m_tree_Mtt->Branch("mHadTop_AfterChi2andKF" , &m_mHadTop_AfterChi2andKF, "mHadTop_AfterChi2andKF/F");
  m_tree_Mtt->Branch("mtt_AfterChi2andKF"     , &m_mtt_AfterChi2andKF    , "mtt_AfterChi2andKF/F");

  // Index of selected particles inside respective collection for mtt computation
  m_tree_Mtt->Branch("selectedLeptonIndex"        , &m_selectedLeptonIndex       , "selectedLeptonIndex/I");
  m_tree_Mtt->Branch("selectedLeptonicBIndex"     , &m_selectedLeptonicBIndex    , "selectedLeptonicBIndex/I");
  m_tree_Mtt->Branch("selectedHadronicBIndex"     , &m_selectedHadronicBIndex    , "selectedHadronicBIndex/I");
  m_tree_Mtt->Branch("selectedHadronicFirstJetIndex"  , &m_selectedHadronicFirstJetIndex  , "selectedHadronicFirstJetIndex/I");
  m_tree_Mtt->Branch("selectedHadronicSecondJetIndex" , &m_selectedHadronicSecondJetIndex , "selectedHadronicSecondJetIndex/I");


  m_tree_Mtt->Branch("trigger_passed", &m_trigger_passed, "trigger_passed/O");

  /// Analysis settings (you define them in your python script)

  // If the setting is not defined the method returns -1. In this
  // case we set a default value for the cut, in order to avoid
  // unwanted crash

  // Main options

  settings->getSetting("doUseBTaginChi2", m_MAIN_doUseBTag);
//    ? m_MAIN_doUseBTag = (static_cast<bool>(settings->getSetting("doUseBTaginChi2")))
//    : m_MAIN_doUseBTag = false;
  settings->getSetting("doChoiceWKF", m_MAIN_doKF);
//    ? m_MAIN_doKF = (static_cast<bool>(settings->getSetting("doChoiceWKF")))
//    : m_MAIN_doKF = false;
  settings->getSetting("doSyst", m_MAIN_doSyst);
//    ? m_MAIN_doSyst = (static_cast<bool>(settings->getSetting("doSyst")))
//    : m_MAIN_doSyst = false;
  settings->getSetting("systvalue", m_MAIN_systvalue);
//    ? m_MAIN_systvalue = settings->getSetting("systvalue")
//    : m_MAIN_systvalue = 0;
  settings->getSetting("doSemiMu", m_MAIN_doSemiMu);
//    ? m_MAIN_doSemiMu = settings->getSetting("doSemiMu")
//    : m_MAIN_doSemiMu = false;

  // VertexSel()
  settings->getSetting("VTX_Ndof_Min", m_VTX_NDof_Min);
//    ? m_VTX_NDof_Min = settings->getSetting("VTX_Ndof_Min") // Value from the joboption
//    : m_VTX_NDof_Min = 0;                                   // Default val


  // METSel()
  settings->getSetting("MET_Pt_Min", m_MET_Pt_Min);
//    ? m_MET_Pt_Min = settings->getSetting("MET_Pt_Min")
//    : m_MET_Pt_Min = 0;


  // MuonSel()
  settings->getSetting("MU_Pt_min_loose", m_MU_Pt_min_loose);
//    ? m_MU_Pt_min_loose = settings->getSetting("MU_Pt_min_loose")
//    : m_MU_Pt_min_loose = 0;
  settings->getSetting("MU_Eta_max_loose", m_MU_Eta_max_loose);
//    ? m_MU_Eta_max_loose = settings->getSetting("MU_Eta_max_loose")
//    : m_MU_Eta_max_loose = 0;
  settings->getSetting("MU_Iso_min", m_MU_Iso_min);
//    ? m_MU_Iso_min = settings->getSetting("MU_Iso_min")
//    : m_MU_Iso_min = 0;
  settings->getSetting("MU_Pt_min", m_MU_Pt_min);
//    ? m_MU_Pt_min = settings->getSetting("MU_Pt_min")
//    : m_MU_Pt_min = 0;
  settings->getSetting("MU_Eta_max", m_MU_Eta_max);
//    ? m_MU_Eta_max = settings->getSetting("MU_Eta_max")
//    : m_MU_Eta_max = 0;
  settings->getSetting("MU_normChi2_max", m_MU_normChi2_max);
//    ? m_MU_normChi2_max = settings->getSetting("MU_normChi2_max")
//    : m_MU_normChi2_max = 0;
  settings->getSetting("MU_nValTrackHits_min", m_MU_nValTrackHits_min);
//    ? m_MU_nValTrackHits_min = settings->getSetting("MU_nValTrackHits_min")
//    : m_MU_nValTrackHits_min = 0;
  settings->getSetting("MU_nMatches_min", m_MU_nMatches_min);
//    ? m_MU_nMatches_min = settings->getSetting("MU_nMatches_min")
//    : m_MU_nMatches_min = 0;
  settings->getSetting("MU_nValPixHits_min", m_MU_nValPixHits_min);
//    ? m_MU_nValPixHits_min = settings->getSetting("MU_nValPixHits_min")
//    : m_MU_nValPixHits_min = 0;
  settings->getSetting("MU_dB_min", m_MU_dB_min);
//    ? m_MU_dB_min = settings->getSetting("MU_dB_min")
//    : m_MU_dB_min = 0;
  settings->getSetting("MU_ePt_min", m_MU_ePt_min);
//    ? m_MU_ePt_min = settings->getSetting("MU_ePt_min")
//    : m_MU_ePt_min = 0;
  settings->getSetting("MU_eEta_max", m_MU_eEta_max);
//    ? m_MU_eEta_max = settings->getSetting("MU_eEta_max")
//    : m_MU_eEta_max = 0;
  settings->getSetting("MU_eEtaW_min", m_MU_eEtaW_min);
//    ? m_MU_eEtaW_min = settings->getSetting("MU_eEtaW_min")
//    : m_MU_eEtaW_min = 0;
  settings->getSetting("MU_eEtaW_max", m_MU_eEtaW_max);
//    ? m_MU_eEtaW_max = settings->getSetting("MU_eEtaW_max")
//    : m_MU_eEtaW_max = 0;
  settings->getSetting("MU_eIso_min", m_MU_eIso_min);
//    ? m_MU_eIso_min = settings->getSetting("MU_eIso_min")
//    : m_MU_eIso_min = 0;


  // ElectronSel()
  settings->getSetting("ELE_Iso_min", m_ELE_Iso_min);
//    ? m_ELE_Iso_min = settings->getSetting("ELE_Iso_min")
//    : m_ELE_Iso_min = 0;
  settings->getSetting("ELE_Pt_min", m_ELE_Pt_min);
//    ? m_ELE_Pt_min = settings->getSetting("ELE_Pt_min")
//    : m_ELE_Pt_min = 0;
  settings->getSetting("ELE_Eta_max", m_ELE_Eta_max);
//    ? m_ELE_Eta_max = settings->getSetting("ELE_Eta_max")
//    : m_ELE_Eta_max = 0;
  settings->getSetting("ELE_Zmass", m_ELE_Zmass);
//    ? m_ELE_Zmass = settings->getSetting("ELE_Zmass")
//    : m_ELE_Zmass = 0;
  settings->getSetting("ELE_Zwin", m_ELE_Zwin);
//    ? m_ELE_Zwin = settings->getSetting("ELE_Zwin")
//    : m_ELE_Zwin = 0;
  settings->getSetting("ELE_dB_min", m_ELE_dB_min);
//    ? m_ELE_dB_min = settings->getSetting("ELE_dB_min")
//    : m_ELE_dB_min = 0;


  // JetSel()
  settings->getSetting("JET_Pt_min", m_JET_Pt_min);
//    ? m_JET_Pt_min = settings->getSetting("JET_Pt_min")
//    : m_JET_Pt_min = 0;
  settings->getSetting("JET_Eta_max", m_JET_Eta_max);
//    ? m_JET_Eta_max = settings->getSetting("JET_Eta_max")
//    : m_JET_Eta_max = 0;

  settings->getSetting("JET_btag_CSVL_min", m_JET_btag_CSVL_min);
//    ? m_JET_btag_CSVL_min = settings->getSetting("JET_btag_CSVL_min")
//    : m_JET_btag_CSVL_min = 0;
  settings->getSetting("JET_btag_CSVM_min", m_JET_btag_CSVM_min);
//    ? m_JET_btag_CSVM_min = settings->getSetting("JET_btag_CSVM_min")
//    : m_JET_btag_CSVM_min = 0;
  settings->getSetting("JET_btag_CSVT_min", m_JET_btag_CSVT_min);
//    ? m_JET_btag_CSVT_min = settings->getSetting("JET_btag_CSVT_min")
//    : m_JET_btag_CSVT_min = 0;
  //(settings->getSetting("JET_btag_TCHPL_min") != -1)
  //  ? m_JET_btag_TCHPL_min = settings->getSetting("JET_btag_TCHPL_min")
  //  : m_JET_btag_TCHPL_min = 0;
  //(settings->getSetting("JET_btag_TCHPM_min") != -1)
  //  ? m_JET_btag_TCHPM_min = settings->getSetting("JET_btag_TCHPM_min")
  //  : m_JET_btag_TCHPM_min = 0;
  settings->getSetting("JET_btag_TCHPT_min", m_JET_btag_TCHPT_min);
//    ? m_JET_btag_TCHPT_min = settings->getSetting("JET_btag_TCHPT_min")
//    : m_JET_btag_TCHPT_min = 0;
  //(settings->getSetting("JET_btag_SSVHEM_min") != -1)
  //  ? m_JET_btag_SSVHEM_min = settings->getSetting("JET_btag_SSVHEM_min")
  //  : m_JET_btag_SSVHEM_min = 0;
  //(settings->getSetting("JET_btag_SSVHPT_min") != -1)
  //  ? m_JET_btag_SSVHPT_min = settings->getSetting("JET_btag_SSVHPT_min")
  //  : m_JET_btag_SSVHPT_min = 0;


  // Triggers
  m_trigger = "";
  settings->getSetting<std::string>("trigger", m_trigger);

  if (m_trigger.length() > 0) {
    m_trigger_regex = boost::regex(m_trigger, boost::regex_constants::optimize);
  }

  TString fname = "kfparams_semilept.dat";

  // Kinfit()
  settings->getSetting("W_mass", m_w);
//    ? m_w = settings->getSetting("W_mass")
//    : m_w = 0;
  settings->getSetting("Top_mass", m_t);
//    ? m_t = settings->getSetting("Top_mass")
//    : m_t = 0;
  settings->getSetting("W_mass_err", m_we);
//    ? m_we = settings->getSetting("W_mass_err")
//    : m_we = 0;
  settings->getSetting("Top_mass_err", m_te);
//    ? m_te = settings->getSetting("Top_mass_err")
//    : m_te = 0;
  settings->getSetting("b_mass", m_b);
//    ? m_b = settings->getSetting("b_mass")
//    : m_b = 0;

  m_KinFit = new KinFit(fname, m_w, m_t, m_b, m_we, m_te);

  //put the complete path to the JEScorr.root file if you run on lyon tier3
  //otherwise remember to put it as an "additional_input_files" in your crab.cfg
  if (m_MAIN_doSyst)
  {
    //TFile* JES_unc_file = new TFile("/gridgroup/cms/Mtt/JEScorr.root");
    //histo_uncertainty= (TH2D*)JES_unc_file->Get("JESfactorJ");
    //jecUnc = new JetCorrectionUncertainty("/gridgroup/cms/Mtt/Jec11_V1_AK5PF_Uncertainty.txt");

  }

}

mtt_analysis_new::~mtt_analysis_new()
{
  delete jecUnc; // Always safe
  delete m_KinFit;
}



//
// First you start with the different physics objects selections
//


// Vertices
int mtt_analysis_new::VertexSel()
{
  int n_vtx = m_vertex->getSize();

  if (!n_vtx)
    return 2;

  for (int i = 0; i < n_vtx; ++i)
  {
    if (m_vertex->getVtxIsFake(i)) continue;
    if (m_vertex->getVtxNdof(i) < m_VTX_NDof_Min) continue;

    return 1;
  }

  return 2;
}


// MET
int mtt_analysis_new::METSel()
{
  m_mtt_MET = m_MET->getMETLorentzVector(0)->Pt();
  if (m_mtt_MET > m_MET_Pt_Min) {
    return 1;
  }

  return 3;
}


int mtt_analysis_new::MuonSel()
{
  int goodmuidx = -1;

  int n_mu = m_muon->getSize();

  if (!n_mu)
    return 4;

  for (int i = 0; i < n_mu; i++)
  {
    if (! m_muon->getMuisGlobal(i))
      continue;

    TLorentzVector *muP = m_muon->getMuLorentzVector(i);

    if (fabs(muP->Pt()) < m_MU_Pt_min_loose)
      continue;
    if (fabs(muP->Eta()) > m_MU_Eta_max_loose)
      continue;


    //add the vertex cut!
    //get the 3-momentum
    //    mu3P = muP->Vect();
    //    Mupass2Dcut=Make2DCut(mu3P,m_jet,MuDRmin,MupTrelmin);
    //    if(Mupass2Dcut==0) continue;
    //for the moment (EPS 2011) isolation cut and not 2D
    float muIso = (m_muon->getMupfChargedHadronIso(i) + m_muon->getMupfNeutralHadronIso(i) + m_muon->getMupfPhotonIso(i)) / fabs(muP->Pt());

    if (muIso > m_MU_Iso_min)
      continue;

    ++m_mtt_NLooseGoodMuons;


    if (m_muon->getMuisTracker(i) != 1)                            continue;
    if (fabs(muP->Pt()) < m_MU_Pt_min)                             continue;
    if (fabs(muP->Eta()) > m_MU_Eta_max)                           continue;
    if (m_muon->getMunormChi2(i) >= m_MU_normChi2_max)             continue;
    if (m_muon->getMunValTrackerHits(i) < m_MU_nValTrackHits_min)  continue;
    if (m_muon->getMunMatches(i) < m_MU_nMatches_min)              continue;
    if (m_muon->getMunValPixelHits(i) < m_MU_nValPixHits_min)      continue;
    if (fabs(m_muon->getMudB(i)) > m_MU_dB_min)                    continue;

    m_mtt_MuonPt[m_mtt_NGoodMuons]   = muP->Pt();
    m_mtt_MuRelIso[m_mtt_NGoodMuons] = muIso;

    ++m_mtt_NGoodMuons;
    goodmuidx = i;
  }

  // No muons
  if (m_mtt_NGoodMuons == 0)
    return 4;

  if (m_mtt_NGoodMuons == 1 && m_mtt_NLooseGoodMuons > 1)
    return 6;

  if (m_mtt_NGoodMuons > 1)
    return 5;


  //electron veto for semimu channel
  int n_ele = m_electron->getSize();

  if (n_ele)
  {
    for (int i = 0; i < n_ele; ++i)
    {
      TLorentzVector *eP = m_electron->getEleLorentzVector(i);

      if (fabs(eP->Pt()) < m_MU_ePt_min)       continue;
      if (fabs(eP->Eta()) > m_MU_eEta_max)     continue;
      if (fabs(eP->Eta()) > m_MU_eEtaW_min
          && fabs(eP->Eta()) < m_MU_eEtaW_max) continue;

      float eleIso = (m_electron->getElepfChargedHadronIso(i) +
          m_electron->getElepfNeutralHadronIso(i) +
          m_electron->getElepfPhotonIso(i)) / eP->Pt();

      if (eleIso > m_MU_eIso_min)
        continue;

      // One good muon, but one good electron too.
      return 7;
    }
  }

  m_refLept = m_muon->getMuLorentzVector(goodmuidx);
  m_selectedLeptonIndex = goodmuidx;

  return 1;
}

int mtt_analysis_new::ElectronSel()
{
  int goodelidx = -1;

  int n_ele = m_electron->getSize();

  if (!n_ele) return 4;

  for (int i = 0; i < n_ele; i++)
  {
    TLorentzVector *eP = m_electron->getEleLorentzVector(i);

    if (fabs(eP->Pt()) < m_ELE_Pt_min)                continue;
    if (fabs(eP->Eta()) > m_ELE_Eta_max)              continue;
    if (fabs(m_electron->getEledB(i)) > m_ELE_dB_min) continue;

    float eleIso = (m_electron->getElepfChargedHadronIso(i) +
        m_electron->getElepfNeutralHadronIso(i) +
        m_electron->getElepfPhotonIso(i)) / eP->Pt();

    if (eleIso > m_ELE_Iso_min)
      continue;

    m_mtt_ElectronPt[m_mtt_NGoodElectrons] = eP->Pt();
    m_mtt_ElRelIso[m_mtt_NGoodElectrons]   = eleIso;

    ((m_electron->getEleHyperTight1MC(i) & 5) == 5)
      ? m_mtt_HyperTight1MC[m_mtt_NGoodElectrons] = 1
      : m_mtt_HyperTight1MC[m_mtt_NGoodElectrons] = 0;

    m_mtt_NGoodElectrons++;

    // Only one good electron per event
    if (m_mtt_NGoodElectrons > 1)
      return 5;

    goodelidx = i;
  }


  // No electrons? bye bye
  if (m_mtt_NGoodElectrons == 0)
    return 4;

  // Z veto
  TLorentzVector *firsteP = m_electron->getEleLorentzVector(goodelidx);

  for (int i = 0; i < n_ele; i++)
  {
    // Don't build invariant mass with ourselves
    if (i == goodelidx)
      continue;

    TLorentzVector *secondeP = m_electron->getEleLorentzVector(i);
    if (fabs(secondeP->Pt()) < m_MU_ePt_min)       continue;
    if (fabs(secondeP->Eta()) > m_MU_eEta_max)     continue;
    if (fabs(secondeP->Eta()) > m_MU_eEtaW_min
        && fabs(secondeP->Eta()) < m_MU_eEtaW_max) continue;

    float eleIso = (m_electron->getElepfChargedHadronIso(i) +
        m_electron->getElepfNeutralHadronIso(i) +
        m_electron->getElepfPhotonIso(i)) / secondeP->Pt();

    if (eleIso > m_ELE_Iso_min)
      continue;

    // Too close to Z mass?
    if (fabs((*firsteP + *secondeP).M() - m_ELE_Zmass) < m_ELE_Zwin)
      return 6;
  }


  //muon veto for semie channel (loose selection)
  int n_mu = m_muon->getSize();

  if (n_mu)
  {
    for (int i = 0; i < n_mu; i++)
    {
      if (!m_muon->getMuisGlobal(i))
        continue;

      TLorentzVector *muP = m_muon->getMuLorentzVector(i);

      if (fabs(muP->Pt()) < m_MU_Pt_min_loose)   continue;
      if (fabs(muP->Eta()) > m_MU_Eta_max_loose) continue;

      float muIso = (m_muon->getMupfChargedHadronIso(i) +
          m_muon->getMupfNeutralHadronIso(i) +
          m_muon->getMupfPhotonIso(i)) / fabs(muP->Pt());

      if (muIso > m_MU_Iso_min)
        continue;

      // We don't want a good muon if we have a good electron
      return 7;
    }
  }


  m_refLept = m_electron->getEleLorentzVector(goodelidx);
  m_selectedLeptonIndex = goodelidx;

  return 1;
}


int mtt_analysis_new::JetSel()
{
  AllJetsPt       = 0.;
  m_selJetsIds.clear();

  int n_jet = m_jet->getSize();

  if (! n_jet)
    return 8;

  for (int i = 0; i < n_jet; i++)
  {
    TLorentzVector *jetP = m_jet->getJetLorentzVector(i);

    if (fabs(jetP->Pt()) < m_JET_Pt_min) continue;
    if (fabs(jetP->Eta()) >  m_JET_Eta_max) continue;

    m_mtt_JetEta[m_mtt_NJets] = jetP->Eta();
    m_mtt_JetPt[m_mtt_NJets]  = jetP->Pt();
    ++m_mtt_NJets;

    AllJetsPt += jetP->Pt();

    if (m_mtt_NJets < 9) // Count the number of btagged jets in the selected jets
    {
      m_selJetsIds.push_back(i);
      if ((m_jet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVL_min)
        ++m_mtt_NBtaggedJets_CSVL;
      if ((m_jet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVM_min)
        ++m_mtt_NBtaggedJets_CSVM;
      if ((m_jet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVT_min)
        ++m_mtt_NBtaggedJets_CSVT;
//      if ((m_jet->getJetBTagProb_TCHP(i)) > m_JET_btag_TCHPL_min)
//        ++m_mtt_NBtaggedJets_TCHPL;
//      if ((m_jet->getJetBTagProb_TCHP(i)) > m_JET_btag_TCHPM_min)
//        ++m_mtt_NBtaggedJets_TCHPM;
      if ((m_jet->getJetBTagProb_TCHP(i)) > m_JET_btag_TCHPT_min)
        ++m_mtt_NBtaggedJets_TCHPT;
//      if ((m_jet->getJetBTagProb_SSVHE(i)) > m_JET_btag_SSVHEM_min)
//        ++m_mtt_NBtaggedJets_SSVHEM;
//      if ((m_jet->getJetBTagProb_SSVHP(i)) > m_JET_btag_SSVHPT_min)
//        ++m_mtt_NBtaggedJets_SSVHPT;
    }

    if (m_mtt_NJets == 1) m_mtt_1stjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
    if (m_mtt_NJets == 2) m_mtt_2ndjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
    if (m_mtt_NJets == 3) m_mtt_3rdjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
    if (m_mtt_NJets == 4) m_mtt_4thjetpt = m_mtt_JetPt[m_mtt_NJets - 1];
  }

  // We need at least 4 good jets
  if (m_mtt_NJets < 4)
    return 8;

  return 1;
}



int mtt_analysis_new::Make2DCut(TVector3 lept3P, float cutDR, float cutPtrel)
{
  pass2Dcut   = 0;
  minjetpt2D  = 30.;
  DrMin       = std::numeric_limits<float>::infinity();
  pTRel       = std::numeric_limits<float>::infinity();
  costheta    = std::numeric_limits<float>::infinity();

  int n_jet = m_jet->getSize();

  if (!n_jet)
    return 0;

  //loop over the jets to calculate variables for 2D cut

  for (int ij = 0; ij < n_jet; ij++)
  {
    TLorentzVector *jetP2D = m_jet->getJetLorentzVector(ij);

    if (fabs(jetP2D->Pt()) < minjetpt2D)
      continue;

    //get the 3-momentum

    jet3P2D = jetP2D->Vect();
    if ((jet3P2D.DeltaR(lept3P)) < DrMin)
    {
      DrMin = jet3P2D.DeltaR(lept3P);
      costheta = ((lept3P.Px() * jetP2D->Px() + lept3P.Py() * jetP2D->Py() + lept3P.Pz() * jetP2D->Pz()) / (lept3P.Mag() * jetP2D->P()));
      pTRel = lept3P.Mag() * sqrt(1. - pow(costheta, 2));
    }
  }

  if (DrMin > cutDR || pTRel > cutPtrel)
    pass2Dcut = 1;

  return pass2Dcut;
}

#define CHECK_RES_AND_RETURN(res) \
  if (res != 1) { \
    m_mtt_isSel = res; \
    return res; \
  }

#define   MTT_TRIGGER_NOT_FOUND   1000

int mtt_analysis_new::mtt_Sel(bool do_MC_, EventExtractor * event, HLTExtractor* HLT, MCExtractor * MC, MuonExtractor *muon, ElectronExtractor *electron, JetExtractor *jet, METExtractor *MET, VertexExtractor *vertex, const edm::EventSetup& iSetup)
{
  reset();

  m_refLept  = nullptr;

  m_vertex   = vertex;
  m_MET      = MET;
  m_muon     = muon;
  m_electron = electron;
  m_jet      = jet;
  m_event    = event;

  if (!m_trigger_regex.empty()) {
    std::vector<std::string>& paths = *HLT->getPaths();

    m_trigger_passed = false;
    for (std::string& path: paths) {
      if (regex_match(path, m_trigger_regex)) {
        //std::cout << "Matched trigger: " << path << std::endl;
        m_trigger_passed = true;
        break;
      }
    }
  } else {
    m_trigger_passed = true;
  }

  if (do_MC_)
  {
    m_MC = MC;
    MCidentification();
  }

  if (m_MAIN_doSyst)
  {

    if (! jecUnc) {
      // Get jet uncertainty from GlobalTag
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      iSetup.get<JetCorrectionsRecord>().get("AK5PFchs",JetCorParColl); 
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      jecUnc = new JetCorrectionUncertainty(JetCorPar);
    }

    SystModifJetsAndMET(1, jecUnc);
  }

  m_nPU = m_event->nPU();

  int res = VertexSel();
  CHECK_RES_AND_RETURN(res);

  res = METSel();
  CHECK_RES_AND_RETURN(res);

  if (m_MAIN_doSemiMu) {
    res  = MuonSel();
  } else {
    res = ElectronSel();
  }
  CHECK_RES_AND_RETURN(res);

  res = JetSel();
  CHECK_RES_AND_RETURN(res);

  m_mtt_isSel = 1;

  // We selected a candidate (and m_refLept)
  loopOverCombinations(do_MC_);

  return 1;
}



void mtt_analysis_new::loopOverCombinations(bool do_MC_)
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

          if (do_MC_)
            m_mtt_OneMatchedCombi = match_MC(c_j1, c_j2, c_j3, c_j4, 0);

          int res = m_KinFit->ReadObjects(*m_jet->getJetLorentzVector(c_j3),
              *m_jet->getJetLorentzVector(c_j4),
              *m_jet->getJetLorentzVector(c_j1),
              *m_refLept,
              *m_MET->getMETLorentzVector(0),
              *m_jet->getJetLorentzVector(c_j2),
              m_MAIN_doSemiMu
              );

          if (!res)
            return; // We will never get anything with this event

          if (m_MAIN_doKF) //use the kinfit to choose the best jet pairing
          {
            (m_KinFit->Fit()) // Do the kinfit converged
              ? fitchi2 = m_KinFit->GetKFChi2()
              : fitchi2 = std::numeric_limits<double>::infinity();
          }
          else  //else use the chi2 to chose the best combination and apply the kinfit only to this latter
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
    m_KinFit->ReadObjects(*m_jet->getJetLorentzVector(bestj3),
        *m_jet->getJetLorentzVector(bestj4),
        *m_jet->getJetLorentzVector(bestj1),
        *m_refLept,
        *m_MET->getMETLorentzVector(0),
        *m_jet->getJetLorentzVector(bestj2),
        m_MAIN_doSemiMu
        );

    const TLorentzVector& measuredLepton = m_KinFit->GetMeasuredLepton();
    const TLorentzVector& measuredNeutrino = m_KinFit->GetMeasuredNeutrino();
    const TLorentzVector& measuredLeptonicB = m_KinFit->GetMeasuredLeptonicBJet();
    const TLorentzVector& measuredHadronicB = m_KinFit->GetMeasuredHadronicBJet();
    const TLorentzVector& measuredHadronicFirstJet = m_KinFit->GetMeasuredFirstLightJet();
    const TLorentzVector& measuredHadronicSecondJet = m_KinFit->GetMeasuredSecondLightJet();

    /**
     * Compute Mtt after doing KinFit
     */
    m_mHadW_AfterChi2   = (measuredHadronicFirstJet + measuredHadronicSecondJet).M();
    m_mLepTop_AfterChi2 = (measuredLepton + measuredNeutrino + measuredLeptonicB).M();
    m_mHadTop_AfterChi2 = (measuredHadronicFirstJet + measuredHadronicSecondJet + measuredHadronicB).M();
    m_mtt_AfterChi2     = (measuredLepton + measuredNeutrino + measuredLeptonicB + measuredHadronicB + measuredHadronicFirstJet + measuredHadronicSecondJet).M();

    /*
       std::cout << "Mt lepton: " << m_mLepTop_AfterChi2 << std::endl;
       std::cout << "Mt hadronic: " << m_mHadTop_AfterChi2 << std::endl;
       std::cout << "Mtt: " << m_mtt_AfterChi2 << std::endl;
       */

    (m_KinFit->Fit()) // Do the kinfit converged
      ? fitchi2 = m_KinFit->GetKFChi2()
      : fitchi2 = std::numeric_limits<double>::infinity();
  }

  if (do_MC_)
    m_mtt_IsBestSolMatched = match_MC(bestj1, bestj2, bestj3, bestj4, 0);

  m_mtt_KFChi2 = fitchi2;

  if (!m_MAIN_doKF)
    m_mtt_BestSolChi2 = minfitchi2;


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

void mtt_analysis_new::MCidentification()
{
  nEle    = 0;
  nMu     = 0;
  nTau    = 0;
  nNuEle  = 0;
  nNuMu   = 0;
  nNuTau  = 0;
  nQuarkb = 0;
  nW      = 0;
  nTop    = 0;
  bool sameW   = false;
  bool sameTop = false;
  idxW.clear();
  idxTop.clear();
  Top.clear();

  int n_MC = m_MC->getSize();

  if (!n_MC)
    return;

  for (int i = 0; i < n_MC ; ++i)
  {
    sameW = false;
    sameTop = false;
    if (fabs(m_MC->getType(m_MC->getMom1Index(m_MC->getMom1Index(i)))) == 6)
    {
      if (fabs(m_MC->getType(m_MC->getMom1Index(i))) == 24)
      {
        /// Count the number of leptons and neutrinos from Top->W
        if (fabs(m_MC->getType(i)) == 11) ++nEle;   //Electron from WTop
        if (fabs(m_MC->getType(i)) == 13) ++nMu;    //Muon	   from WTop
        if (fabs(m_MC->getType(i)) == 15) ++nTau;   //Tau	   from WTop
        if (fabs(m_MC->getType(i)) == 12) ++nNuEle; //NuEle    from WTop
        if (fabs(m_MC->getType(i)) == 14) ++nNuMu;  //NuMu	   from WTop
        if (fabs(m_MC->getType(i)) == 16) ++nNuTau; //NuTau    from WTop

        /// Count the number of W (with no double counting) from Top
        for (unsigned int j = 0; j < idxW.size(); j++)
        {
          if (idxW.size() != 0 && idxW[j] == m_MC->getMom1Index(i))
          {
            sameW = true;
            break;
          }
        }
        idxW.push_back(m_MC->getMom1Index(i));
        if (!sameW)
        {
          nW++;
        }

        /// Get Top info (with no double counting)
        for (unsigned int k = 0; k < idxTop.size(); k++)
        {
          if (idxTop.size() != 0 && idxTop[k] == m_MC->getMom1Index(m_MC->getMom1Index(i)))
          {
            sameTop = true;
            break;
          }
        }
        idxTop.push_back(m_MC->getMom1Index(m_MC->getMom1Index(i)));
        if (!sameTop)
        {
          TLorentzVector TL_Top(m_MC->getPx(m_MC->getMom1Index(m_MC->getMom1Index(i))),
              m_MC->getPy(m_MC->getMom1Index(m_MC->getMom1Index(i))),
              m_MC->getPz(m_MC->getMom1Index(m_MC->getMom1Index(i))),
              m_MC->getE(m_MC->getMom1Index(m_MC->getMom1Index(i))));


          Top.push_back(TL_Top);
          nTop++;
        }
      }
    }

    /// Count the number of b quark from Top
    if (fabs(m_MC->getType(i)) == 5 && fabs(m_MC->getType(m_MC->getMom1Index(i))) == 6)
    {
      nQuarkb++; //Quark b from Top
    }
  }

  if (nEle == 1 && nNuEle == 1 && nMu == 0 && nNuMu == 0 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nW == 2 && nTop == 2)
  {
    m_MC_channel = 1;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 1 && nNuMu == 1 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nW == 2 && nTop == 2)
  {
    m_MC_channel = 2;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 0 && nNuMu == 0 && nTau == 1 && nNuTau == 1 && nQuarkb > 1 && nW == 2 && nTop == 2)
  {
    m_MC_channel = 3;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 0 && nNuMu == 0 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nW == 2 && nTop == 2)
  {
    m_MC_channel = 4;
  }

  if (nEle == 2 && nNuEle == 2 && nMu == 0 && nNuMu == 0 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nW == 2 && nTop == 2)
  {
    m_MC_channel = 5;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 2 && nNuMu == 2 && nTau == 0 && nNuTau == 0 && nQuarkb > 1 && nW == 2 && nTop == 2)
  {
    m_MC_channel = 6;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 0 && nNuMu == 0 && nTau == 2 && nNuTau == 2 && nQuarkb > 1 && nW == 2 && nTop == 2)
  {
    m_MC_channel = 7;
  }

  if (nEle == 1 && nNuEle == 1 && nMu == 1 && nNuMu == 1 && nTau == 0 && nNuTau == 0 && nQuarkb == 2 && nW == 2 && nTop == 2)
  {
    m_MC_channel = 8 ;
  }

  if (nEle == 1 && nNuEle == 1 && nMu == 0 && nNuMu == 0 && nTau == 1 && nNuTau == 1 && nQuarkb == 2 && nW == 2 && nTop == 2)
  {
    m_MC_channel = 9 ;
  }

  if (nEle == 0 && nNuEle == 0 && nMu == 1 && nNuMu == 1 && nTau == 1 && nNuTau == 1 && nQuarkb == 2 && nW == 2 && nTop == 2)
  {
    m_MC_channel = 10;
  }

  if (Top.size() >= 2) {
    m_MC_mtt = (Top[0] + Top[1]).M();
  } else {
    m_MC_mtt = -1;
  }

  // Extract index of semi-leptonic event, and store them in tree. Useful if you want to know how many jets you have selected right
  if (false) {
    std::cout << "New event" << std::endl;
    for (int i = 0; i < n_MC; i++) {
      std::cout << "\tType: " << m_MC->getType(i) << std::endl;
    }
  }

  bool keepEvent = true;
  for (int i = 0; i < n_MC; i++) {
    if (m_MC->getMom1Index(i) == -1 || m_MC->getMom1Index(m_MC->getMom1Index(i)) == -1)
      continue;

    // Look only event coming (directly / indirectly) from a top
    if (fabs(m_MC->getType(m_MC->getMom1Index(i))) == ID_T ||
        fabs(m_MC->getType(m_MC->getMom1Index(m_MC->getMom1Index(i)))) == ID_T)  {

      int type = (int) fabs(m_MC->getType(i));
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

          m_leptonicWIndex = m_MC->getMom1Index(i);
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

  if (! keepEvent)
    return;

  // Reorder B jet indexes
  if (m_MC->getMom1Index(m_leptonicBIndex) != m_leptonicWIndex) {
    // Wrong combinaison, swap
    std::swap(m_leptonicBIndex, m_hadronicBIndex);
  }
}

int mtt_analysis_new::match_MC(int idxJetbH, int idxJetbL, int idxJet1,	int idxJet2,
    int idxLepton)
{
  if (
      /// Ask if Jet b hadronique  come from a b and top
      fabs(m_MC->getType(m_jet->getJetMCIndex(idxJetbH))) == 5 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_jet->getJetMCIndex(idxJetbH)))) == 6 &&

      /// Ask if Jet b leptonique  come from a b and top
      fabs(m_MC->getType(m_jet->getJetMCIndex(idxJetbL))) == 5 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_jet->getJetMCIndex(idxJetbL)))) == 6 &&

      /// Ask if jet 1,2 come from light quark and W and top
      fabs(m_MC->getType(m_jet->getJetMCIndex(idxJet1))) < 5 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_jet->getJetMCIndex(idxJet1)))) == 24 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_MC->getMom1Index(m_jet->getJetMCIndex(idxJet1))))) == 6 &&
      fabs(m_MC->getType(m_jet->getJetMCIndex(idxJet2))) < 5 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_jet->getJetMCIndex(idxJet2)))) == 24 &&
      fabs(m_MC->getType(m_MC->getMom1Index(m_MC->getMom1Index(m_jet->getJetMCIndex(idxJet2))))) == 6
     )
  {
    return 1;
  }
  else
  {
    return 0;
  }
}


void mtt_analysis_new::fillTree()
{
  m_tree_Mtt->Fill();
}

void mtt_analysis_new::SystModifJetsAndMET(int SystType, JetCorrectionUncertainty *jecUnc)
{
  double unc = 0.;
  double corr = 0.;
  double met_corr_x = 0.;
  double met_corr_y = 0.;
  unsigned int nJets = m_jet->getSize();
  TLorentzVector *myMET = m_MET->getMETLorentzVector(0);

  if (SystType == 1)   // JES
  {
    for (unsigned int iJet = 0; iJet < nJets; iJet++)
    {
      TLorentzVector *myJet = m_jet->getJetLorentzVector(iJet);

      jecUnc->setJetEta(myJet->Eta());
      jecUnc->setJetPt(myJet->Pt()); // here you must use the CORRECTED jet pt

      unc = (SystType == 1) ? jecUnc->getUncertainty(true) : jecUnc->getUncertainty(false);
      corr = fabs(unc);

      //use corrected jet pt for met correction
      met_corr_x += myJet->Px() * (m_MAIN_systvalue * corr);
      met_corr_y += myJet->Py() * (m_MAIN_systvalue * corr);

      m_jet->setJetLorentzVector(iJet, myJet->E() * (1. + m_MAIN_systvalue * corr), myJet->Px() * (1. + m_MAIN_systvalue * corr), myJet->Py() * (1. + m_MAIN_systvalue * corr), myJet->Pz() * (1. + m_MAIN_systvalue * corr));
    }

    m_MET->setMETLorentzVector(0, myMET->E(), myMET->Px() - met_corr_x, myMET->Py() - met_corr_y, myMET->Pz());
  }
}


// Here we just reset the ROOTtree parameters

void mtt_analysis_new::reset()
{
  jecUnc = NULL;

  m_mtt_isSel = 0;
  m_mtt_IsBestSolMatched = -1;
  m_mtt_OneMatchedCombi = 0;
  m_mtt_BestSolChi2 = -1.;
  m_mtt_KFChi2 = -1.;
  m_mtt_AfterChi2andKF = -1.;
  m_mtt_NumComb = 0;
  m_mLepTop_AfterChi2andKF = -1.;
  m_mHadTop_AfterChi2andKF = -1.;

  m_mtt_AfterChi2     = -1.;
  m_mLepTop_AfterChi2 = -1.;
  m_mHadTop_AfterChi2 = -1.;
  m_mHadW_AfterChi2   = -1.;

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
  m_leptonicWIndex = -1;

  m_firstJetIndex = -1;
  m_secondJetIndex = -1;
  
  m_trigger_passed = false;
}
