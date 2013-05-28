#pragma once

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
#include <boost/regex.hpp>

#include <Extractors/PatExtractor/interface/ExtractorPlugin.h>

#include <TTree.h>
#include <TLorentzVector.h>

class AnalysisSettings;
class EventExtractor;
class MuonExtractor;
class ElectronExtractor;
class METExtractor;
class VertexExtractor;
class KinFit;
class HLTExtractor;
class PatExtractor;
class MCExtractor;
class JetMETExtractor;

namespace edm {
  class EventSetup;
}

namespace patextractor {

  class mtt_analysis: public Plugin {
    public:
      mtt_analysis(const edm::ParameterSet& cmsswSettings);
      ~mtt_analysis();

      // TTbar selection
      virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, PatExtractor& extractor);
      virtual void analyze(const edm::EventSetup& iSetup, PatExtractor& extractor);

    private:

      //Selections
      int MuonSel();
      int ElectronSel();
      int JetSel();
      int VertexSel();
      int METSel();

      void loopOverCombinations();

      //makes the 2D cut for a given lepton
      //int Make2DCut(TVector3 lept3P, float cutDR, float cutPtrel);

      /// MC event channel identification
      int match_MC(int idxJetbH, int idxJetbL, int idxJet1, int idxJet2, int idxLepton);

      void MCidentification();
      void reset();
      void fillTree();

      bool isBJet(unsigned int index); 

      bool   m_MAIN_doUseBTag;
      //bool   m_MAIN_doKF;
      bool   m_MAIN_doSemiMu;

      TTree*  m_tree_Mtt;
      KinFit* m_KinFit;

      std::shared_ptr<EventExtractor> m_event;
      std::shared_ptr<MCExtractor>    m_MC;

      std::shared_ptr<VertexExtractor> m_vertex;

      //std::shared_ptr<METExtractor>   m_MET;
      float m_MET_Pt_Min;

      std::shared_ptr<MuonExtractor> m_muon;
      float m_MU_Iso_max;
      float m_MU_Pt_min;
      float m_MU_Eta_max;

      std::shared_ptr<MuonExtractor> m_muon_loose;
      float m_MU_Pt_min_loose;
      float m_MU_Eta_max_loose;
      float m_MU_Iso_max_loose;

      std::shared_ptr<ElectronExtractor> m_electron;
      float m_ELE_Pt_min;
      float m_ELE_Eta_max;
      float m_ELE_Iso_max;

      std::shared_ptr<ElectronExtractor> m_electron_loose;
      float m_ELE_Pt_min_loose;
      float m_ELE_Eta_max_loose;
      float m_ELE_Iso_max_loose;

      std::shared_ptr<JetMETExtractor> m_jetMet;
      float m_JET_Pt_min;
      float m_JET_Eta_max;
      //float m_jet_btag_tchel_min;
      //float m_jet_btag_tchem_min;
      //float m_jet_btag_tchet_min;
      //float m_JET_btag_TCHPL_min;
      //float m_JET_btag_TCHPM_min;
      float m_JET_btag_TCHPT;
      //float m_JET_btag_SSVHEM_min;
      //float m_JET_btag_SSVHPT_min;
      float m_JET_btag_CSVL;
      float m_JET_btag_CSVM;
      float m_JET_btag_CSVT;

      TLorentzVector* m_refLept;

      // Triggers
      bool m_trigger_passed;

      //MC stuff
      int m_MC_channel;
      float m_MC_mtt;
      int m_nPU; // number of interactions

      // Indexes of gen particle in a semi-lept event
      int m_leptonIndex;
      int m_neutrinoIndex;

      int m_leptonicBIndex;
      int m_hadronicBIndex;
      int m_leptonicTopIndex;

      int m_firstJetIndex;
      int m_secondJetIndex;

      float m_MC_hadronicWMass;
      float m_MC_leptonicWMass;
      float m_MC_hadronicTopMass;
      float m_MC_leptonicTopMass;
      float m_MC_pt_tt;

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
      //float m_mtt_KFChi2;

      /*float m_mtt_AfterChi2andKF;
        float m_mLepTop_AfterChi2andKF;
        float m_mHadTop_AfterChi2andKF;*/

      float m_mtt_AfterChi2;
      float m_pt_tt_AfterChi2;
      float m_mLepTop_AfterChi2;
      float m_mHadTop_AfterChi2;
      float m_mHadW_AfterChi2;
      float m_mLepW_AfterChi2;

      float m_lepTopPt_AfterChi2;
      float m_lepTopEta_AfterChi2;
      float m_hadTopPt_AfterChi2;
      float m_hadTopEta_AfterChi2;

      // Indexes of selected particles for mtt computation
      int m_selectedLeptonIndex;
      int m_selectedLeptonicBIndex;
      int m_selectedHadronicBIndex;
      int m_selectedHadronicFirstJetIndex;
      int m_selectedHadronicSecondJetIndex;



      int m_mtt_NGoodMuons;
      int m_mtt_NLooseGoodMuons;
      float m_mtt_MuonPt[20];
      float m_mtt_MuonEta[20];
      float m_mtt_MuRelIso[20];
      float m_mtt_2DDrMin[20];
      float m_mtt_2DpTrel[20];

      int m_mtt_NGoodElectrons;
      float m_mtt_ElectronPt[20];
      float m_mtt_ElectronEta[20];
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

      float m_b_tagging_efficiency;

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


      // Kin Fit
      /// maximal number of iterations to be performed for the fit
      unsigned int maxNrIter_;
      /// maximal chi2 equivalent
      double maxDeltaS_;
      /// maximal deviation for contstraints
      double maxF_;
      unsigned int jetParam_;
      unsigned int lepParam_;
      unsigned int metParam_;
      /// constrains
      std::vector<unsigned> constraints_;
      double mW_;
      double mTop_;
      /// scale factors for jet energy resolution
      std::vector<double> jetEnergyResolutionScaleFactors_;
      std::vector<double> jetEnergyResolutionEtaBinning_;
      /// config-file-based object resolutions
      std::vector<edm::ParameterSet> udscResolutions_;
      std::vector<edm::ParameterSet> bResolutions_;
      std::vector<edm::ParameterSet> lepResolutions_;
      std::vector<edm::ParameterSet> metResolutions_;

      float m_weight;
      float m_weight_error_low;
      float m_weight_error_high;

      bool m_is_neutrino_pz_corrected;

      // Cut ; -1 event drop before arriving to this cut ; 0 cut failed, 1 cut passed
      int m_pass_vertex_cut;
      int m_pass_met_cut;
      int m_pass_lepton_cut;
      int m_pass_jet_cut;
  };

}
