#pragma once

#include "../interface/SuperBaseExtractor.h"

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/Candidate/interface/CandidateFwd.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include <iostream>


/**
 * \class MCExtractor
 * \brief Stores information about generator-level particles
 * 
 * 
 */
class MCExtractor: public SuperBaseExtractor
{
public:
    /**
     * \brief Constructor to be used when running the "extractor + analyses" chain
     * 
     * Configures the extractor to store information in a tree.
     */
    MCExtractor(const std::string& name, const edm::ParameterSet& settings);
    
    /**
     * \brief Constructor to be used when only analyses are run
     * 
     * Cofigures the extractor to read information stored in a tree created by a previous run of
     * this extractor.
     */
    MCExtractor(const std::string& name, const edm::ParameterSet& settings, TFile* srcFile);
    
public:
    /// Declares which data the extractor is going to read from an EDM event
    virtual void doConsumes(edm::ConsumesCollector&& collector) override;
    
    /**
     * \brief Processes current event in the source EDM file
     * 
     * The last parameter is not used.
     */
    virtual void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup,
     MCExtractor*) override;
    
    /// Initialization to be performed for each event
    void reset();
    
    /// Fills main tree with properties of particles
    void fillTree();
    
    /// Reads a new event from the source file
    void getInfo(int ievt);
    
    int getSize() const;
    int getStatus(int index);
    int getType(int index);
    float getPx(int index);
    float getPy(int index);
    float getPz(int index);
    float getE(int index);
    float getMom1Index(int index);
    TLorentzVector* getP4(int index);
    int getPatIndex(int index) const;
    
private:
    void constructGeneration(int gene, int npart);
    
private:
    TTree* m_tree_MC;
    edm::InputTag m_genParticleTag;
    edm::EDGetTokenT<edm::View<reco::GenParticle>> m_genParticleToken;
    bool m_doJpsi;
    
    static const Int_t m_MCs_MAX = 1000;
    
    Int_t m_n_MCs;
    
    TClonesArray *m_MC_lorentzvector;
    Int_t   m_MC_status[m_MCs_MAX];
    Int_t   m_MC_type[m_MCs_MAX];
    Int_t   m_MC_imot1[m_MCs_MAX];
    Int_t   m_MC_generation[m_MCs_MAX];
    Float_t	m_MC_vx[m_MCs_MAX];
    Float_t	m_MC_vy[m_MCs_MAX];
    Float_t	m_MC_vz[m_MCs_MAX];
        
    // Arrays below contain redundant information. They should be removed in future
    Int_t   m_MC_index[m_MCs_MAX];
    Int_t   m_MC_imot2[m_MCs_MAX];
    Float_t m_MC_E[m_MCs_MAX];
    Float_t m_MC_px[m_MCs_MAX];
    Float_t m_MC_py[m_MCs_MAX];
    Float_t m_MC_pz[m_MCs_MAX];
    Float_t	m_MC_eta[m_MCs_MAX];
    Float_t	m_MC_phi[m_MCs_MAX];
    Bool_t m_MC_LeptonFromTop[m_MCs_MAX];
    Bool_t m_MC_LeptonFromAntiTop[m_MCs_MAX];
    Bool_t m_MC_JPsiFromTop[m_MCs_MAX];
    Bool_t m_MC_JPsiFromAntiTop[m_MCs_MAX];
};
