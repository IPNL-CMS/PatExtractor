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

    JetMETExtractor(const std::string& name, const std::string& met_name, const edm::ParameterSet& config);
    JetMETExtractor(const std::string& name, const std::string& met_name, TFile *a_file);
    virtual ~JetMETExtractor();

    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* m_MC); 

    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Jet& part, int index); 
    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::MET& part, int index); 

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
    void setMETLorentzVector(int idx, float E, float Px, float Py, float Pz);

    // Jet ID
    bool isPFJetLoose(const pat::Jet& jet);

    void correctMETWithTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets);
    void correctJets(pat::JetCollection& jets, const edm::Event& iEvent, const edm::EventSetup& iSetup);

    double getResCorrFactor(const pat::Jet& jet);
    void correctJetsMETresolution(pat::JetCollection& jets, pat::MET& met);

    double getSysShifCorrFactorX(const int Nvtx);
    double getSysShifCorrFactorY(const int Nvtx);
    void correctMETWithSysShift(const edm::Event& event, pat::MET& met);

    void doJESSystematics(pat::JetCollection& jets, pat::MET& met);

    void extractRawJets(pat::JetCollection& jets);

  private:

    TTree* m_tree_jet;

    bool mCorrectJets;
    bool mCorrectSysShiftMet;
    std::string mJetCorrectorLabel;
    GreaterByPt<pat::Jet> mSorter;
    bool mRedoTypeI;
    bool mDoJER;
    int  mJERSign;
    int  mJESSign;

    static const int 	m_jets_MAX       = 200;

    TClonesArray* m_jet_lorentzvector;
    float	m_jet_vx[m_jets_MAX];
    float	m_jet_vy[m_jets_MAX];
    float	m_jet_vz[m_jets_MAX];
    int	m_jet_chmult[m_jets_MAX];
    float	m_jet_chmuEfrac[m_jets_MAX];
    float	m_jet_chemEfrac[m_jets_MAX];
    float	m_jet_chhadEfrac[m_jets_MAX];
    float	m_jet_nemEfrac[m_jets_MAX];
    float	m_jet_nhadEfrac[m_jets_MAX];

    //float	m_jet_btag_BjetProb[m_jets_MAX];
    //float	m_jet_btag_SSVHE[m_jets_MAX];
    //float	m_jet_btag_SSVHP[m_jets_MAX];
    //float	m_jet_btag_TCHE[m_jets_MAX];
    float	m_jet_btag_jetProb[m_jets_MAX];
    float	m_jet_btag_TCHP[m_jets_MAX];
    float	m_jet_btag_CSV[m_jets_MAX];

    int  m_jet_MCIndex[m_jets_MAX];

    // MET

    edm::InputTag m_metTag;
    TTree* m_tree_met;
    TClonesArray* m_met_lorentzvector;

    ScaleFactorCollection m_scaleFactors;

    std::shared_ptr<JetCorrectionUncertainty> jecUncertainty;
    std::string mJecFilename;
};

#endif 
