#ifndef KVFEXTRACTOR_H
#define KVFEXTRACTOR_H

/**
 * KVFExtractor
 * \brief: Base class for extracting jet info
 */


//Include PAT info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <CommonTools/Utils/interface/PtComparator.h>

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/View.h"
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>
#include "CMGTools/External/interface/PileupJetIdentifier.h"
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "../interface/BaseExtractor.h"
#include "../interface/MCExtractor.h"
#include <Extractors/PatExtractor/interface/ScaleFactor.h>
#include <Extractors/PatExtractor/interface/ScaleFactorService.h>

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "../interface/MuScleFitCorrector.h"

class JetCorrectionUncertainty;

class KVFExtractor: public BaseExtractor<pat::Jet>
{

  public:

    KVFExtractor(const std::string& name_jpsi, const std::string& name_d0, const std::string& name_jet, std::shared_ptr<ScaleFactorService> sf, const edm::ParameterSet& config);
    KVFExtractor(const std::string& name_jpsi, const std::string& name_d0, const std::string& name_jet, std::shared_ptr<ScaleFactorService> sf, TFile *a_file);
    virtual ~KVFExtractor();
    virtual void beginJob();

    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* m_MC); 

    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Jet& part, int index) {}

    void getInfo(int ievt); 

    void reset();
    void fillTree(); 

    /* 
     * Need to be there !
     */
    
    virtual const reco::Candidate* getGenParticle(const pat::Jet& jet) {
      return jet.genParton();
    }

    virtual void setGenParticleIndex(int genParticleIndex, int index) {
    }

    virtual float getMCMatchDeltaR() {
      return 0.4;
    }

    virtual float getMCMatchDPtRel() {
      return 3.0;
    }

    virtual std::vector<int> getPdgIds() {
      return {1, 2, 3, 4, 5, 21};
    }

    virtual TLorentzVector getP4(const pat::Jet& object) {
      TLorentzVector p4;
      p4.SetPxPyPzE(object.px(), object.py(), object.pz(), object.energy());

      return p4;
    }

    virtual double getLesUp(double eta) {
      // from AN2014_170_v2, table 8, DATA-MC column
      double lesUnc = 1.;
      if (fabs(eta) < 1.5)
        lesUnc += 0.05/100.;
      else if (fabs(eta) < 2.4) 
        lesUnc += 0.10/100.;
      return lesUnc;
    }

    virtual double getLesDown(double eta) {
      // from AN2014_170_v2, table 8, DATA-MC column
      double lesUnc = 1.;
      if (fabs(eta) < 1.5)
        lesUnc -= 0.05/100.;
      else if (fabs(eta) < 2.4) 
        lesUnc -= 0.10/100.;
      return lesUnc;
    }


    // Jet ID
    bool isPFJetLoose(const pat::Jet& jet);
    
    void correctJets(pat::JetCollection& jets, const edm::Event& iEvent, const edm::EventSetup& iSetup);

    double getResCorrFactor(const pat::Jet& jet);
    void correctJetsResolution(pat::JetCollection& jets);

    double getSysShifCorrFactorX(const int Nvtx);
    double getSysShifCorrFactorY(const int Nvtx);

    void doJESSystematics(pat::JetCollection& jets);

    void extractRawJets(pat::JetCollection& jets);


  private:

    double m_muJpsiMinPt;
    double m_jpsiMassMin;
    double m_jpsiMassMax;
    int m_nTrD0Max;
    double m_trSumMinPt;
    double m_trUnfoldMinPt;

    bool mCorrectJets;
    bool mUseGlobalTagForJEC;
    std::string mJecPayload;
    std::string mJecJetAlgo;
    std::string mJetCorrectorLabel;
    GreaterByPt<pat::Jet> mSorterJets;
    bool mDoJER;
    bool mDoLooseJetID;
    int  mJERSign;
    int  mJESSign;
    
    FactorizedJetCorrector* mTxtCorrector; 

    std::shared_ptr<JetCorrectionUncertainty> jecUncertainty;
    std::string mJecFilename;

    // Jpsi

    edm::InputTag m_vertexTag;

    TTree* m_tree_jpsi;

    static const int 	m_jpsi_MAX  = 10;

    std::string fitParametersFile;
    MuScleFitCorrector* corrector;
    
    int           m_jpsi_size;
    int           m_jpsi_indjet[m_jpsi_MAX];
    float         m_jpsi_jet_btag_CSV[m_jpsi_MAX];
    TClonesArray* m_jpsi_jet_lorentzvector;
    ScaleFactorCollection m_jpsi_jet_scaleFactors;
    int           m_jpsi_indpf1[m_jpsi_MAX];
    int           m_jpsi_indpf2[m_jpsi_MAX];
    TClonesArray* m_jpsipf_lorentzvector;
    TClonesArray* m_jpsipf_lesup_lorentzvector;
    TClonesArray* m_jpsipf_lesdown_lorentzvector;
    TClonesArray* m_jpsipf_raw_lorentzvector;
    TClonesArray* m_jpsipf_mu1_lorentzvector;
    TClonesArray* m_jpsipf_mu1_lesup_lorentzvector;
    TClonesArray* m_jpsipf_mu1_lesdown_lorentzvector;
    TClonesArray* m_jpsipf_mu1_raw_lorentzvector;
    TClonesArray* m_jpsipf_mu2_lorentzvector;
    TClonesArray* m_jpsipf_mu2_lesup_lorentzvector;
    TClonesArray* m_jpsipf_mu2_lesdown_lorentzvector;
    TClonesArray* m_jpsipf_mu2_raw_lorentzvector;
    TClonesArray* m_jpsikvf_lorentzvector;
    TClonesArray* m_jpsikvf_lesup_lorentzvector;
    TClonesArray* m_jpsikvf_lesdown_lorentzvector;
    TClonesArray* m_jpsikvf_raw_lorentzvector;
    TClonesArray* m_jpsikvf_mu1_lorentzvector;
    TClonesArray* m_jpsikvf_mu1_lesup_lorentzvector;
    TClonesArray* m_jpsikvf_mu1_lesdown_lorentzvector;
    TClonesArray* m_jpsikvf_mu1_raw_lorentzvector;
    int           m_jpsikvf_mu1_pdgid[m_jpsi_MAX];
    std::map<std::string, ScaleFactorCollection> m_jpsikvf_mu1_muon_scaleFactors;
    TClonesArray* m_jpsikvf_mu2_lorentzvector;
    TClonesArray* m_jpsikvf_mu2_lesup_lorentzvector;
    TClonesArray* m_jpsikvf_mu2_lesdown_lorentzvector;
    TClonesArray* m_jpsikvf_mu2_raw_lorentzvector;
    int           m_jpsikvf_mu2_pdgid[m_jpsi_MAX];
    std::map<std::string, ScaleFactorCollection> m_jpsikvf_mu2_muon_scaleFactors;
    float         m_jpsikvf_vx[m_jpsi_MAX];
    float         m_jpsikvf_vy[m_jpsi_MAX];
    float         m_jpsikvf_vz[m_jpsi_MAX];
    bool          m_jpsikvf_vtxvalid[m_jpsi_MAX];
    float         m_jpsikvf_vtxchi2[m_jpsi_MAX];
    float         m_jpsikvf_ndf[m_jpsi_MAX];
    float         m_jpsikvf_L3D[m_jpsi_MAX];
    float         m_jpsikvf_sigmaL3D[m_jpsi_MAX];
    float         m_jpsikvf_L3DoverSigmaL3D[m_jpsi_MAX];

    // D0

    GreaterByPt<reco::PFCandidate> mSorterPFs;

    TTree* m_tree_mujet;

    static const int 	m_mujet_MAX = 20;
    static const int 	m_d0_MAX  = 100;
    static const int 	m_D0_MAX  = m_d0_MAX*m_mujet_MAX;
    static const int 	m_tr_MAX  = 10;
    static const int 	m_Tr_MAX  = m_tr_MAX*m_mujet_MAX;
    static const int 	m_mu_MAX  = 3;
    static const int 	m_Mu_MAX  = m_mu_MAX*m_mujet_MAX;
    static const int 	m_unfold_tr_MAX = 500;
    static const int  m_unfold_Tr_MAX = m_unfold_tr_MAX*m_mujet_MAX;
    
    int           m_mujet_size;
    int           m_mujet_mu_size;
    int           m_mujet_tr_size;
    int           m_mujet_unfold_tr_size;
    int           m_mujet_d0_size;
    float         m_mujet_jet_btag_CSV[m_mujet_MAX];
    TClonesArray* m_mujet_jet_lorentzvector;
    ScaleFactorCollection m_mujet_jet_scaleFactors;
    int           m_mujet_mu[m_mujet_MAX];
    int           m_mujet_mu_indmujet[m_Mu_MAX]; 
    TClonesArray* m_mujet_mu_lorentzvector; 
    int           m_mujet_mu_pdgid[m_Mu_MAX]; 
    float         m_mujet_mu_iso[m_Mu_MAX]; 
    int           m_mujet_ntr[m_mujet_MAX];
    float         m_mujet_sump[m_mujet_MAX];
    float         m_mujet_sumpt[m_mujet_MAX];
    float         m_mujet_sumvecp[m_mujet_MAX];
    int           m_mujet_tr_indmujet[m_Tr_MAX]; 
    TClonesArray* m_mujet_tr_lorentzvector; 
    int           m_mujet_tr_pdgid[m_Tr_MAX]; 
    int           m_mujet_nd0[m_mujet_MAX];
    int           m_mujet_d0kvf_indmujet[m_D0_MAX]; 
    TClonesArray* m_mujet_d0pf_lorentzvector;
    TClonesArray* m_mujet_d0kvf_lorentzvector;
    TClonesArray* m_mujet_d0kvf_pion_lorentzvector;
    int           m_mujet_d0kvf_pion_pdgid[m_D0_MAX];
    TClonesArray* m_mujet_d0kvf_kaon_lorentzvector;
    int           m_mujet_d0kvf_kaon_pdgid[m_D0_MAX];
    float         m_mujet_d0kvf_vx[m_D0_MAX];
    float         m_mujet_d0kvf_vy[m_D0_MAX];
    float         m_mujet_d0kvf_vz[m_D0_MAX];
    bool          m_mujet_d0kvf_vtxvalid[m_D0_MAX];
    float         m_mujet_d0kvf_vtxchi2[m_D0_MAX];
    float         m_mujet_d0kvf_ndf[m_D0_MAX];
    float         m_mujet_d0kvf_L3D[m_D0_MAX];
    float         m_mujet_d0kvf_sigmaL3D[m_D0_MAX];
    float         m_mujet_d0kvf_L3DoverSigmaL3D[m_D0_MAX];

    int           m_mujet_unfold_indmujet[m_unfold_Tr_MAX]; 
    float         m_mujet_unfold_tr_recopt[m_unfold_Tr_MAX];
    float         m_mujet_unfold_tr_recoeta[m_unfold_Tr_MAX];
    float         m_mujet_unfold_tr_genpt[m_unfold_Tr_MAX];
    float         m_mujet_unfold_tr_geneta[m_unfold_Tr_MAX];
    float         m_mujet_unfold_tr_dr[m_unfold_Tr_MAX];
    float         m_mujet_unfold_mu_recopt[m_unfold_Tr_MAX];
    float         m_mujet_unfold_mu_recoeta[m_unfold_Tr_MAX];
    float         m_mujet_unfold_mu_genpt[m_unfold_Tr_MAX];
    float         m_mujet_unfold_mu_geneta[m_unfold_Tr_MAX];
    float         m_mujet_unfold_mu_dr[m_unfold_Tr_MAX];
    
};

#endif 
