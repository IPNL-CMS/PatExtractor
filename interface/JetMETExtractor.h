#ifndef JETEXTRACTOR_H
#define JETEXTRACTOR_H

/**
 * JetMETExtractor
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

class JetMETExtractor: public BaseExtractor<pat::Jet>
{

  public:

    JetMETExtractor(const std::string& name, const edm::ParameterSet& config);
    JetMETExtractor(const std::string& name, const edm::ParameterSet& config, TFile *a_file);
    virtual ~JetMETExtractor();
    virtual void doConsumes(edm::ConsumesCollector&& collector);
    virtual void beginJob(bool isInAnalysisMode);

    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* m_MC); 

    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Jet& part, int index, const pat::JetRef& ref); 
    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Jet& part, int index) {}
    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::MET& part, int index); 
    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const reco::PFCandidate& part, int index);

    void getInfo(int ievt); 

    void reset();
    void fillTree(); 

    virtual const reco::Candidate* getGenParticle(const pat::Jet& jet) {
      return jet.genParton();
    }

    /** Set the index of the gen particle matched with this jet.
     * The index is relative to MCExtractor array
     */
    virtual void setGenParticleIndex(int genParticleIndex, int index) {
      m_jet_MCIndex[index] = genParticleIndex;
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

    TLorentzVector *getP4(int jetidx) {return (TLorentzVector*)m_jet_lorentzvector->At(jetidx);}
    TLorentzVector *getGenP4(int jetidx) {return (TLorentzVector*)m_genjet_lorentzvector->At(jetidx);}
    TLorentzVector *getRawP4(int jetidx) {return (TLorentzVector*)m_rawjet_lorentzvector->At(jetidx);}


    float getJetBTagProb_CSV(int jetidx) const { return m_jet_btag_CSV[jetidx]; }
    float getJetBTagProb_TCHP(int jetidx) const { return m_jet_btag_TCHP[jetidx]; }

    //float getJetBTagProb_SSVHP(int jetidx) {return m_jet_btag_SSVHP[jetidx];}
    //float getJetBTagProb_TCHE(int jetidx) {return m_jet_btag_TCHE[jetidx];}

    int getJetMCIndex(int jetidx) {return m_jet_MCIndex[jetidx];}

    /**
     * This method is need for inline JEC. We need to be able to change jet pt on-the-fly
     */
    void setJetLorentzVector(int index, TLorentzVector* vector)
    {
      if (index >= m_jet_lorentzvector->GetEntries())
        return;

      new((*m_jet_lorentzvector)[index]) TLorentzVector(*vector);
    }

    void setJetLorentzVector(int index, double energy, double px, double py, double pz)
    {
      TLorentzVector* v = new TLorentzVector(px, py, pz, energy);
      setJetLorentzVector(index, v);
      delete v;
    }

    ScaleFactor getScaleFactor(int index) const {
      return m_scaleFactors.at(index);
    }

    TLorentzVector *getMETLorentzVector(int metidx) {return (TLorentzVector*)m_met_lorentzvector->At(metidx);}
    TLorentzVector *getGenMETLorentzVector(int metidx) {return (TLorentzVector*)m_genmet_lorentzvector->At(metidx);}
    void setMETLorentzVector(int idx, float E, float Px, float Py, float Pz);

    int getNumberOfUnclusteredParticles() {return m_unclustered_particle_lorentzvector->GetEntriesFast();}
    
    TLorentzVector *getUnclusteredParticleLorentzVector(int unclustpcidx) {return (TLorentzVector*)m_unclustered_particle_lorentzvector->At(unclustpcidx);}

    // Jet ID
    bool isPFJetLoose(const pat::Jet& jet);
    int isPFJetLoose(int muidx) const { return m_jet_isPFJetLoose[muidx]; }
    
    // PU jet ID
    int getPuJetFullId(int muidx) const { return m_jet_puJetFullId[muidx]; }
    int getPuJetCutBasedId(int muidx) const { return m_jet_puJetCutBasedId[muidx]; }


    void correctMETWithTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets);
    void correctJets(pat::JetCollection& jets, const edm::Event& iEvent, const edm::EventSetup& iSetup);

    double getResCorrFactor(const pat::Jet& jet);
    void correctJetsMETresolution(pat::JetCollection& jets, pat::MET& met);

    double getSysShifCorrFactorX(const int Nvtx);
    double getSysShifCorrFactorY(const int Nvtx);
    void correctMETWithSysShift(const edm::Event& event, pat::MET& met);

    void doJESSystematics(pat::JetCollection& jets, pat::MET& met);

    void extractRawJets(pat::JetCollection& jets);

    int getAlgoPartonFlavor(int index) const {
      return m_jet_algo_parton_flavor[index];
    }

    int getPhysicsPartonFlavor(int index) const {
      return m_jet_physics_parton_pdgid[index];
    }

  private:

    TTree* m_tree_jet;
    edm::InputTag m_rawMetTag;
    edm::InputTag m_particleFlowTag;
    edm::InputTag m_primaryVerticesTag;
    edm::InputTag m_rhoTag;

    edm::EDGetTokenT<pat::METCollection> m_metToken;
    edm::EDGetTokenT<pat::METCollection> m_rawMetToken;
    edm::EDGetTokenT<reco::PFCandidateCollection> m_particleFlowToken;
    edm::EDGetTokenT<reco::VertexCollection> m_primaryVerticesToken;
    edm::EDGetTokenT<double> m_rhoToken;

    bool mCorrectJets;
    bool mUseGlobalTagForJEC;
    bool mSaveUnclusteredParticles;
    std::string mJecPayload;
    std::string mJecJetAlgo;
    bool mCorrectSysShiftMet;
    std::string mJetCorrectorLabel;
    GreaterByPt<pat::Jet> mSorter;
    bool mRedoTypeI;
    bool mDoJER;
    bool mDoLooseJetID;
    int  mJERSign;
    int  mJESSign;
    
    FactorizedJetCorrector* mTxtCorrector; 

    static const int 	m_jets_MAX       = 200;

    TClonesArray* m_jet_lorentzvector;
    TClonesArray* m_genjet_lorentzvector;
    TClonesArray* m_rawjet_lorentzvector;

    float	m_jet_vx[m_jets_MAX];
    float	m_jet_vy[m_jets_MAX];
    float	m_jet_vz[m_jets_MAX];
    int   m_jet_chmult[m_jets_MAX];
    float	m_jet_chmuEfrac[m_jets_MAX];
    float	m_jet_chemEfrac[m_jets_MAX];
    float	m_jet_chhadEfrac[m_jets_MAX];
    float	m_jet_nemEfrac[m_jets_MAX];
    float	m_jet_nhadEfrac[m_jets_MAX];
    int   m_jet_isPFJetLoose[m_jets_MAX];

    //float	m_jet_btag_BjetProb[m_jets_MAX];
    //float	m_jet_btag_SSVHE[m_jets_MAX];
    //float	m_jet_btag_SSVHP[m_jets_MAX];
    //float	m_jet_btag_TCHE[m_jets_MAX];
    float	m_jet_btag_jetProb[m_jets_MAX];
    float	m_jet_btag_TCHP[m_jets_MAX];
    float	m_jet_btag_CSV[m_jets_MAX];
    
    //PuJetId
    int    m_jet_puJetFullId[m_jets_MAX];
    int    m_jet_puJetCutBasedId[m_jets_MAX];
    // is loose wp  : 4
    // is medium wp : 6
    // is tight wp  : 7


    int  m_jet_MCIndex[m_jets_MAX];

    int m_jet_algo_parton_flavor[m_jets_MAX];
    int m_jet_physics_parton_pdgid[m_jets_MAX];

    // MET

    edm::InputTag m_metTag;
    TTree* m_tree_met;
    TClonesArray* m_met_lorentzvector;
    float m_met_sumEt;
    TClonesArray* m_unclustered_particle_lorentzvector;
    TClonesArray* m_genmet_lorentzvector;

    ScaleFactorCollection m_scaleFactors;

    std::shared_ptr<JetCorrectionUncertainty> jecUncertainty;
    std::string mJecFilename;
};

#endif 
