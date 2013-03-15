#ifndef ELECTRONEXTRACTOR_H
#define ELECTRONEXTRACTOR_H

/**
 * ElectronExtractor
 * \brief: Base class for extracting electron info
 */


//Include PAT info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
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

class ElectronExtractor: public BaseExtractor<pat::Electron>
{

  public:

    ElectronExtractor(const std::string& name, const edm::InputTag& tag, bool doTree);
    ElectronExtractor(const std::string& name, TFile* f);
    virtual ~ElectronExtractor();

    virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Electron& object, int index);

    void reset();
    void fillTree(); 
    void getInfo(int ievt); 

    virtual const reco::Candidate* getGenParticle(const pat::Electron& electron) {
      return electron.genLepton();
    }

    virtual void setGenParticleIndex(int genParticleIndex, int index) {
      m_ele_MCIndex[index] = genParticleIndex;
    }

    virtual float getMCMatchDeltaR() {
      return 0.5;
    }

    virtual float getMCMatchDPtRel() {
      return 0.5;
    }

    virtual std::vector<int> getPdgIds() {
      return {11};
    }

    virtual TLorentzVector getP4(const pat::Electron& object) {
      TLorentzVector p4;
      p4.SetPxPyPzE(object.px(), object.py(), object.pz(), object.energy());

      return p4;
    }

    TLorentzVector *getEleLorentzVector(int eidx) {return (TLorentzVector*)m_ele_lorentzvector->At(eidx);}
    float getEledB(int eidx) {return m_ele_dB[eidx];}
    float getElepfChargedHadronIso(int eidx) {return m_ele_pfChargedHadronIso[eidx];}
    float getElepfNeutralHadronIso(int eidx) {return m_ele_pfNeutralHadronIso[eidx];}
    float getElepfPhotonIso(int eidx) {return m_ele_pfPhotonIso[eidx];}

    int getEleMCIndex(int eidx){return m_ele_MCIndex[eidx];}
    int getEleCharge(int eidx) {return m_ele_charge[eidx];}

    /*
    int getEleHyperTight1MC(int eidx){return m_ele_eidHyperTight1MC[eidx];}
    int getEleLooseMC(int eidx)      {return m_ele_eidLooseMC[eidx];}
    int getEleMediumMC(int eidx)     {return m_ele_eidMediumMC[eidx];}
    int getEleSuperTightMC(int eidx) {return m_ele_eidSuperTightMC[eidx];}
    int getEleTightMC(int eidx)      {return m_ele_eidTightMC[eidx];}
    int getEleVeryLooseMC(int eidx)  {return m_ele_eidVeryLooseMC[eidx];}
    int getElepf_evspi(int eidx)     {return m_ele_eidpf_evspi[eidx];}
    int getElepf_evsmu(int eidx)     {return m_ele_eidpf_evsmu[eidx];}
    */

    bool passVetoId(int index) const {
      return m_ele_passVetoID[index];
    }

    bool passLooseId(int index) const {
      return m_ele_passLooseID[index];
    }

    bool passMediumId(int index) const {
      return m_ele_passMediumID[index];
    }

    bool passTightId(int index) const {
      return m_ele_passTightID[index];
    }

    float getRhoCorrectedRelativeIsolation(int index) const {
      return m_ele_rhoCorrectedRelIsolation[index];
    }

    float getSuperClusterEta(int index) const {
      return m_ele_SCEta[index];
    }

    ScaleFactor getScaleFactor(int index) const {
      return m_scaleFactors.at(index);
    }

  private:

    TTree* m_tree_electron;

    static const int 	m_electrons_MAX  = 100;

    float m_deltaR_cut;

    TClonesArray* m_ele_lorentzvector;
    float	m_ele_vx[m_electrons_MAX];
    float	m_ele_vy[m_electrons_MAX];
    float	m_ele_vz[m_electrons_MAX];
    int	m_ele_charge[m_electrons_MAX];
    float m_ele_SCEta[m_electrons_MAX];
    bool m_ele_passConversionVeto[m_electrons_MAX];

    // electron id's

    /// old one
    int   m_ele_eidLoose[m_electrons_MAX]; 
    int   m_ele_eidRobustHighEnergy[m_electrons_MAX]; 
    int   m_ele_eidRobustLoose[m_electrons_MAX]; 
    int   m_ele_eidRobustTight[m_electrons_MAX]; 
    int   m_ele_eidTight[m_electrons_MAX]; 
    int   m_ele_eidpf_evspi[m_electrons_MAX]; 
    int   m_ele_eidpf_evsmu[m_electrons_MAX];

    /// new one
    int m_ele_eidHyperTight1MC[m_electrons_MAX];
    int m_ele_eidLooseMC[m_electrons_MAX];
    int m_ele_eidMediumMC[m_electrons_MAX];
    int m_ele_eidSuperTightMC[m_electrons_MAX];
    int m_ele_eidTightMC[m_electrons_MAX];
    int m_ele_eidVeryLooseMC[m_electrons_MAX];

    /// 2012
    float m_ele_eidMVATrigV0[m_electrons_MAX];
    bool m_ele_passVetoID[m_electrons_MAX];
    bool m_ele_passLooseID[m_electrons_MAX];
    bool m_ele_passMediumID[m_electrons_MAX];
    bool m_ele_passTightID[m_electrons_MAX];
    float m_ele_effectiveArea[m_electrons_MAX];
    float m_ele_rhoCorrectedRelIsolation[m_electrons_MAX]; // Isolation corrected with effective area
    float m_ele_deltaBetaCorrectedRelIsolation[m_electrons_MAX]; // Isolation corrected with delta beta corrections
    float m_ele_relIsolation[m_electrons_MAX]; // Isolation

    //
    int eidBit;
    bool pass ;
    float m_ele_dB[m_electrons_MAX];
    float m_ele_trackIso[m_electrons_MAX];
    float m_ele_ecalIso[m_electrons_MAX];
    float m_ele_hcalIso[m_electrons_MAX];
    float m_ele_pfParticleIso[m_electrons_MAX]; // isolation calculated with all the PFCandidates
    float m_ele_pfChargedHadronIso[m_electrons_MAX]; // isolation calculated with only the charged hadron PFCandidates
    float m_ele_pfNeutralHadronIso[m_electrons_MAX]; // isolation calculated with only the neutral hadron PFCandidates
    float m_ele_pfPhotonIso[m_electrons_MAX]; // Returns the isolation calculated with only the photon PFCandidates
    int   m_ele_numberOfMissedInnerLayer[m_electrons_MAX]; // Access the hit pattern counting (in the Tracker) the number of expected crossed layers  before the first trajectory's hit
    int   m_ele_MCIndex[m_electrons_MAX];

    ScaleFactorCollection m_scaleFactors;

};

#endif 
