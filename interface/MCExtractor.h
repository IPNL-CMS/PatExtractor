#pragma once

#include "../interface/SuperBaseExtractor.h"

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>


/**
 * \class MCExtractor
 * \brief Stores information about generator-level particles
 * 
 * Stores properties of generator-level particles from the final state of the hard interaction
 * before modelling of FSR, as well as prompt electrons, muons, and photons. If configured
 * specially, stores also J/psi decaying to a pair of muons and these muons.
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
    
    /// Destructor
    ~MCExtractor();
    
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
    
    /// Fills main tree with properties of particles
    void fillTree();
    
    /// Reads a new event from the source file
    void getInfo(int ievt);
    
    /// Returns the number of stored particles
    int getSize() const;
    
    /// Returns PDG ID of stored particle with the given index
    int getType(int index) const;
    
    /**
     * \brief Return status code of stored particle with the given index
     * 
     * You should avoid using explicit status codes in an analysis since they have no strict
     * physical meaning.
     */
    int getStatus(int index) const;
    
    /// Returns four-momentum of stored particle with the given index
    TLorentzVector const &getP4(int index) const;
    
    /**
     * \brief Return index of the first stored ancestor of the particle with the given index
     * 
     * The returned index can be used to access the ancestor directly, e.g. obtain its momentum with
     * getP4(getMom1Index(index)). If no ancestors of the particle have been stored, the method
     * returns -1.
     */
    int getMom1Index(int index) const;
    
    /**
     * \depricated Index returned by the getMom1Index method can be used to access ancestor directly
     */
    int getPatIndex(int index) const __attribute__ ((deprecated));
    
    /**
     * \brief Returns x component of momentum of stored particle with the given index
     * 
     * \depricated Access momentum with the getP4 method instead.
     */
    float getPx(int index) const __attribute__ ((deprecated));
    
    /**
     * \brief Returns y component of momentum of stored particle with the given index
     * 
     * \depricated Access momentum with the getP4 method instead.
     */
    float getPy(int index) const __attribute__ ((deprecated));
    
    /**
     * \brief Returns z component of momentum of stored particle with the given index
     * 
     * \depricated Access momentum with the getP4 method instead.
     */
    float getPz(int index) const __attribute__ ((deprecated));
    
    /**
     * \brief Returns energy of stored particle with the given index
     * 
     * \depricated Access energy with the getP4 method instead.
     */
    float getE(int index) const __attribute__ ((deprecated));
    
private:
    /// Name of collection of generator-level particles
    edm::InputTag m_genParticleTag;
    
    /// Token to read the collection of generator-level particles
    edm::EDGetTokenT<edm::View<reco::GenParticle>> m_genParticleToken;
    
    /// Flag controlling if the extractor stores information related to J/psi
    bool m_doJpsi;
    
    /**
     * \brief Output tree
     * 
     * Created by the extractor but not owned by it.
     */
    TTree* m_tree_MC;
    
    /// Maximal number of particles that can be stored
    static const Int_t m_MCs_MAX = 100;
    
    /// Number of particles in the current event
    Int_t m_n_MCs;
    
    /// PDG ID of stored particles
    Int_t m_MC_type[m_MCs_MAX];
    
    /// Status codes of stored particles
    Int_t m_MC_status[m_MCs_MAX];
    
    /**
     * \brief Four-momenta of stored particles
     * 
     * The object is owned by the extractor.
     */
    TClonesArray *m_MC_lorentzvector;
    
    /// Index of stored ancestors of stored particles
    Int_t m_MC_imot1[m_MCs_MAX];
    
    /// Coordinates of vertices of stored particles
    Float_t	m_MC_vx[m_MCs_MAX], m_MC_vy[m_MCs_MAX],	m_MC_vz[m_MCs_MAX];
    
        
    // Arrays below contain redundant information. They should be removed in future
    Int_t m_MC_index[m_MCs_MAX];
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
