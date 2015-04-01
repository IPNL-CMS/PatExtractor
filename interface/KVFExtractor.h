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

#include "../interface/BaseExtractor.h"
#include "../interface/MCExtractor.h"
#include <Extractors/PatExtractor/interface/ScaleFactor.h>

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class JetCorrectionUncertainty;

class KVFExtractor: public BaseExtractor<pat::Jet>
{

  public:

    KVFExtractor(const std::string& name, std::shared_ptr<ScaleFactorService> sf, const edm::ParameterSet& config);
    KVFExtractor(const std::string& name, std::shared_ptr<ScaleFactorService> sf, TFile *a_file);
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

    TTree* m_tree_jpsi;

    bool mCorrectJets;
    bool mUseGlobalTagForJEC;
    std::string mJecPayload;
    std::string mJecJetAlgo;
    std::string mJetCorrectorLabel;
    GreaterByPt<pat::Jet> mSorter;
    bool mDoJER;
    bool mDoLooseJetID;
    int  mJERSign;
    int  mJESSign;
    
    FactorizedJetCorrector* mTxtCorrector; 

    std::shared_ptr<JetCorrectionUncertainty> jecUncertainty;
    std::string mJecFilename;

    // Jpsi

    static const int 	m_jpsi_MAX  = 10;
    
    int           m_jpsi_size;
    int           m_jpsi_indjet[m_jpsi_MAX];
    int           m_jpsi_indpf1[m_jpsi_MAX];
    int           m_jpsi_indpf2[m_jpsi_MAX];
    TClonesArray* m_jpsi_jet_lorentzvector;
    TClonesArray* m_jpsipf_lorentzvector;
    TClonesArray* m_jpsikvf_lorentzvector;
    TClonesArray* m_jpsikvf_mu1_lorentzvector;
    TClonesArray* m_jpsikvf_mu2_lorentzvector;
    float         m_jpsikvf_vx[m_jpsi_MAX];
    float         m_jpsikvf_vy[m_jpsi_MAX];
    float         m_jpsikvf_vz[m_jpsi_MAX];
    bool          m_jpsikvf_vtxvalid[m_jpsi_MAX];
    float         m_jpsikvf_vtxchi2[m_jpsi_MAX];
    float         m_jpsikvf_ndf[m_jpsi_MAX];
    float         m_jpsikvf_L3D[m_jpsi_MAX];
    float         m_jpsikvf_sigmaL3D[m_jpsi_MAX];
    float         m_jpsikvf_L3DoverSigmaL3D[m_jpsi_MAX];
    
};

#endif 
