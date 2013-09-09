#ifndef PFPARTEXTRACTOR_H
#define PFPARTEXTRACTOR_H

/**
 * PFpartExtractor
 * \brief: Base class for extracting PF particles info
 */


//Include PAT info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "../interface/BaseExtractor.h"
#include "../interface/MCExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class PFpartExtractor: public BaseExtractor<reco::PFCandidate>
{

  public:

    PFpartExtractor(const std::string& name, const edm::InputTag& tag, bool doTree);
    PFpartExtractor(const std::string& name, TFile* a_file);
    virtual ~PFpartExtractor();
    
    // Dummy pure virtual function in BaseExtractor
    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const reco::PFCandidate& part, int index);
    // Function doing the real job
    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor);

    virtual const reco::Candidate* getGenParticle(const reco::PFCandidate& part) {
      return NULL;
    }

    virtual void setGenParticleIndex(int genParticleIndex, int index) {
    }

    virtual float getMCMatchDeltaR() {
      return 0.;
    }

    virtual float getMCMatchDPtRel() {
      return 0.;
    }

    virtual std::vector<int> getPdgIds() {
      return {};
    }

    virtual TLorentzVector getP4(const reco::PFCandidate& object) {
      TLorentzVector p4;
      p4.SetPxPyPzE(object.px(), object.py(), object.pz(), object.energy());

      return p4;
    }

    void reset();
    void fillTree();
    void getInfo(int ievt); 

    // Setters/Getters

    TLorentzVector *getPFLorentzVector(int muidx) {return (TLorentzVector*)m_pf_lorentzvector->At(muidx);}

  private:

    TTree* m_tree_pfpart;

    static const int 	m_pfpart_MAX  = 100;

    int         m_pf_size;
    TClonesArray* m_pf_lorentzvector;
    float	m_pf_vx[m_pfpart_MAX];
    float	m_pf_vy[m_pfpart_MAX];
    float	m_pf_vz[m_pfpart_MAX];
    int	        m_pf_charge[m_pfpart_MAX];
    int         m_pf_pdgid[m_pfpart_MAX];
    int         m_pf_trkLayer[m_pfpart_MAX];
    int         m_pf_pixLayer[m_pfpart_MAX];
    float       m_pf_trknormChi2[m_pfpart_MAX];
    int         m_pf_numberOfChambers[m_pfpart_MAX];
    int         m_pf_numberOfMatchedStations[m_pfpart_MAX];
    bool        m_pf_isGlobalMuon[m_pfpart_MAX];
    bool        m_pf_isTrackerMuon[m_pfpart_MAX];
    bool        m_pf_isStandAloneMuon[m_pfpart_MAX];
    bool        m_pf_isCaloMuon[m_pfpart_MAX];
    bool        m_pf_isPFMuon[m_pfpart_MAX];
    bool        m_pf_isRPCMuon[m_pfpart_MAX];
    
    static const int 	m_jpsi_MAX  = 10;
    
    int           m_jpsi_size;
    int           m_jpsi_indpf1[m_jpsi_MAX];
    int           m_jpsi_indpf2[m_jpsi_MAX];
    TClonesArray* m_jpsimu_lorentzvector;
    TClonesArray* m_jpsiraw_lorentzvector;
    TClonesArray* m_jpsiraw_mu1_lorentzvector;
    TClonesArray* m_jpsiraw_mu2_lorentzvector;
    float         m_jpsiraw_vx[m_jpsi_MAX];
    float         m_jpsiraw_vy[m_jpsi_MAX];
    float         m_jpsiraw_vz[m_jpsi_MAX];
    bool          m_jpsiraw_vtxvalid[m_jpsi_MAX];
    float         m_jpsiraw_vtxchi2[m_jpsi_MAX];
    float         m_jpsiraw_ndf[m_jpsi_MAX];
    float         m_jpsiraw_L3D[m_jpsi_MAX];
    float         m_jpsiraw_sigmaL3D[m_jpsi_MAX];
    float         m_jpsiraw_L3DoverSigmaL3D[m_jpsi_MAX];
    
};

#endif 
