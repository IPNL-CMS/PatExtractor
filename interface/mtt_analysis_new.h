#ifndef MTT_ANALYSIS_NEW_H
#define MTT_ANALYSIS_NEW_H

/*****************************
 * chanel :
 *  0  -> no ttbar decays
 *  1  -> isSemiElectronic
 *  2  -> isSemiMuonic
 *  3  -> isSemiTauic
 *  4  -> isFullHadronic
 *  5  -> isDiElectronic
 *  6  -> isDiMuonic
 *  7  -> isDiTauic
 *  8  -> isElectroMuonic
 *  9  -> isElectroTauic
 *  10 -> isMuoTauic
 ******************************/

#include <vector>

class AnalysisSettings;
class EventExtractor;
class MuonExtractor;
class ElectronExtractor;
class JetExtractor;
class METExtractor;
class VertexExtractor;
class KinFit;

class JetCorrectionUncertainty;

namespace edm {
  class EventSetup;
}

class mtt_analysis_new {
public:
  mtt_analysis_new(AnalysisSettings* settings);
  ~mtt_analysis_new();

  //Selections
  int MuonSel();
  int ElectronSel();
  int JetSel();
  int VertexSel();
  int METSel();

  // TTbar selection
  int mtt_Sel(bool do_MC_, EventExtractor * event, MCExtractor * MC, MuonExtractor *muon, ElectronExtractor *electron, JetExtractor *jet, METExtractor *MET, VertexExtractor *vertex, const edm::EventSetup& iSetup);

  void   loopOverCombinations(bool do_MC_);

  //makes the 2D cut for a given lepton
  int   Make2DCut(TVector3 lept3P, float cutDR, float cutPtrel);


  /// MC event channel identification

  int match_MC(int idxJetbH,
               int idxJetbL,
               int idxJet1,
               int idxJet2,
               int idxLepton);

  void MCidentification();
  void reset();
  void fillTree();

  void SystModifJetsAndMET(int SystType, JetCorrectionUncertainty* jecUnc);

  JetCorrectionUncertainty* jecUnc;

  bool isBJet(unsigned int index) {
    // Use recommanded WP from https://indico.cern.ch/getFile.py/access?contribId=4&resId=2&materialId=slides&confId=195042
    return m_jet->getJetBTagProb_CSV(index) > m_JET_btag_CSVM_min;
  }

private:

  bool   m_MAIN_doUseBTag;
  bool   m_MAIN_doKF;
  bool   m_MAIN_doSyst;
  double m_MAIN_systvalue;
  bool   m_MAIN_doSemiMu;

  TTree*  m_tree_Mtt;
  KinFit* m_KinFit;
  float  m_w;
  float  m_b;
  float  m_t;
  float  m_we;
  float  m_te;

  EventExtractor* m_event;
  MCExtractor*    m_MC;

  VertexExtractor* m_vertex;
  float m_VTX_NDof_Min;

  METExtractor*   m_MET;
  float          m_MET_Pt_Min;

  MuonExtractor* m_muon;
  float m_MU_Pt_min_loose;
  float m_MU_Eta_max_loose;
  float m_MU_Iso_min;
  float m_MU_Pt_min;
  float m_MU_Eta_max;
  float m_MU_normChi2_max;
  float m_MU_nValTrackHits_min;
  float m_MU_nMatches_min;
  float m_MU_nValPixHits_min;
  float m_MU_dB_min;
  float m_MU_ePt_min;
  float m_MU_eEta_max;
  float m_MU_eEtaW_min;
  float m_MU_eEtaW_max;
  float m_MU_eIso_min;

  ElectronExtractor* m_electron;
  float m_ELE_Pt_min;
  float m_ELE_Eta_max;
  float m_ELE_Iso_min;
  float m_ELE_Zmass;
  float m_ELE_Zwin;
  float m_ELE_dB_min;

  JetExtractor* m_jet;
  float m_JET_Pt_min;
  float m_JET_Eta_max;
  //float m_jet_btag_tchel_min;
  //float m_jet_btag_tchem_min;
  //float m_jet_btag_tchet_min;
  //float m_JET_btag_TCHPL_min;
  //float m_JET_btag_TCHPM_min;
  float m_JET_btag_TCHPT_min;
  //float m_JET_btag_SSVHEM_min;
  //float m_JET_btag_SSVHPT_min;
  float m_JET_btag_CSVL_min;
  float m_JET_btag_CSVM_min;
  float m_JET_btag_CSVT_min;

  TLorentzVector* m_refLept;


  //MC stuff


  int m_MC_channel;
  float m_MC_mtt;
  int m_nPU; // number of interactions
  /// Number of lepton/neutrino from Top->W and quark b from Top
  int nEle;
  int nMu;
  int nTau;
  int nNuEle;
  int nNuMu;
  int nNuTau;
  int nQuarkb;
  int nW;
  int nTop;
  std::vector<int> idxW;
  std::vector<int> idxTop;
  std::vector<TLorentzVector> Top;

  std::vector<int> m_selJetsIds;
  float AllJetsPt;

  int SelLeptIdx;
  //Reco stuff



  int m_mtt_isSel;
  int m_mtt_IsBestSolMatched;
  int m_mtt_OneMatchedCombi;

  int m_mtt_NumComb;
  float m_mtt_SolChi2[1000];
  float m_mtt_BestSolChi2;
  float m_mtt_KFChi2;

  float m_mtt_AfterChi2andKF;
  float m_mLepTop_AfterChi2andKF;
  float m_mHadTop_AfterChi2andKF;

  float m_mtt_AfterChi2;
  float m_mLepTop_AfterChi2;
  float m_mHadTop_AfterChi2;



  int m_mtt_NGoodMuons;
  int m_mtt_NLooseGoodMuons;
  float m_mtt_MuonPt[20];
  float m_mtt_MuRelIso[20];
  float m_mtt_2DDrMin[20];
  float m_mtt_2DpTrel[20];

  int m_mtt_NGoodElectrons;
  float m_mtt_ElectronPt[20];
  float m_mtt_ElRelIso[20];
  int m_mtt_HyperTight1MC[20];

  float m_mtt_1stjetpt;
  float m_mtt_2ndjetpt;
  float m_mtt_3rdjetpt;
  float m_mtt_4thjetpt;
  float m_mtt_MET;

  int m_mtt_NJets;
  int m_mtt_NGoodJets;
  //int m_mtt_NBtaggedJets_TCHEL;
  //int m_mtt_NBtaggedJets_TCHEM;
  //int m_mtt_NBtaggedJets_TCHET;
  //int m_mtt_NBtaggedJets_TCHPL;
  //int m_mtt_NBtaggedJets_TCHPM;
  int m_mtt_NBtaggedJets_TCHPT;
  //int m_mtt_NBtaggedJets_SSVHEM;
  //int m_mtt_NBtaggedJets_SSVHPT;
  int m_mtt_NBtaggedJets_CSVL;
  int m_mtt_NBtaggedJets_CSVM;
  int m_mtt_NBtaggedJets_CSVT;

  float m_mtt_GoodJetEta[1000];
  float m_mtt_JetEta[1000];
  float m_mtt_JetPt[1000];

  //generic variables for 2D cut
  int pass2Dcut;

  TVector3 jet3P2D;
  float minjetpt2D;
  float DrMin;
  float pTRel;
  float costheta;

  int nGoodElectrons_veto;
  int Elepass2Dcut_veto;

  //variables for semie selection


  bool itsaZ;
  TVector3 el3P;

  int nGoodMuons_veto;
  //variables for jet selection
  int  nGoodJets;


};

#endif
