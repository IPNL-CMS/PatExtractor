#include "../interface/MCExtractor.h"

#include <algorithm>


using namespace std;


MCExtractor::MCExtractor(const std::string& name, const edm::ParameterSet& settings):
    SuperBaseExtractor(name, settings),
    m_genParticleTag(settings.getParameter<edm::InputTag>("input")),
    m_doJpsi(settings.getParameter<bool>("do_jpsi"))
{
    // Create the output tree. Apparently, at this point it is not associated yet with the output
    //file
    m_tree_MC = new TTree(name.c_str(), "PAT MC info");  
    m_tree_MC->SetAutoSave(0);
    
    
    m_MC_lorentzvector = new TClonesArray("TLorentzVector");
    
    
    // Set up branches
    m_tree_MC->Branch("n_MCs", &m_n_MCs, "n_MCs/I");
    m_tree_MC->Branch("MC_type", &m_MC_type, "MC_type[n_MCs]/I");
    m_tree_MC->Branch("MC_isLastCopy", &m_MC_isLastCopy, "MC_isLastCopy[n_MCs]/O");
    m_tree_MC->Branch("MC_4vector", "TClonesArray", &m_MC_lorentzvector, 1000, 0);
    m_tree_MC->Branch("MC_mot1", &m_MC_imot1, "MC_mot1[n_MCs]/I");
    m_tree_MC->Branch("MC_vx", &m_MC_vx, "MC_vx[n_MCs]/F");
    m_tree_MC->Branch("MC_vy", &m_MC_vy, "MC_vy[n_MCs]/F");
    m_tree_MC->Branch("MC_vz", &m_MC_vz, "MC_vz[n_MCs]/F");
    
    
    // Below are redundant branches, which should be dropped eventually
    m_tree_MC->Branch("MC_index",   &m_MC_index,    "MC_index[n_MCs]/I");
    m_tree_MC->Branch("MC_e",   &m_MC_E,    "MC_e[n_MCs]/F");
    m_tree_MC->Branch("MC_px",  &m_MC_px,   "MC_px[n_MCs]/F");
    m_tree_MC->Branch("MC_py",  &m_MC_py,   "MC_py[n_MCs]/F");
    m_tree_MC->Branch("MC_pz",  &m_MC_pz,   "MC_pz[n_MCs]/F");
    m_tree_MC->Branch("MC_eta", &m_MC_eta,  "MC_eta[n_MCs]/F");
    m_tree_MC->Branch("MC_phi", &m_MC_phi,  "MC_phi[n_MCs]/F");
    
    if (m_doJpsi)
    {
        m_tree_MC->Branch("MC_LeptonFromTop", &m_MC_LeptonFromTop, "m_MC_LeptonFromTop[n_MCs]/O");  
        m_tree_MC->Branch("MC_LeptonFromAntiTop", &m_MC_LeptonFromAntiTop,
         "m_MC_LeptonFromAntiTop[n_MCs]/O");
        m_tree_MC->Branch("MC_JPsiFromTop", &m_MC_JPsiFromTop, "m_MC_JPsiFromTop[n_MCs]/O");  
        m_tree_MC->Branch("MC_JPsiFromAntiTop", &m_MC_JPsiFromAntiTop,
         "m_MC_JPsiFromAntiTop[n_MCs]/O");  
    }
    
    
    m_OK = true;
}


MCExtractor::MCExtractor(const std::string& name, const edm::ParameterSet& settings,
 TFile *srcFile):
    SuperBaseExtractor(name, settings, srcFile),
    m_doJpsi(settings.getParameter<bool>("do_jpsi"))
{
    std::cout << "MCExtractor object is retrieved" << std::endl;
    
    // Read the tree from the source file
    m_OK = false;
    
    m_tree_MC = dynamic_cast<TTree*>(srcFile->Get(name.c_str()));
    
    if (not m_tree_MC)
    {
        std::cerr << "Tree \"" << name << "\" does not exist." << std::endl;
        return;
    }
    
    m_OK = true;
    
    
    m_MC_lorentzvector = new TClonesArray("TLorentzVector");
    
    
    // Set up mandatory branches
    m_tree_MC->SetBranchAddress("n_MCs", &m_n_MCs);
    m_tree_MC->SetBranchAddress("MC_type", &m_MC_type);
    m_tree_MC->SetBranchAddress("MC_4vector", &m_MC_lorentzvector);
    
    
    // Set up optional branches
    TBranch *b;
    
    if ((b = m_tree_MC->FindBranch("MC_status")))
        b->SetAddress(&m_MC_status);
    
    if ((b = m_tree_MC->FindBranch("MC_isLastCopy")))
        b->SetAddress(&m_MC_isLastCopy);
    
    if ((b = m_tree_MC->FindBranch("MC_mot1")))
        b->SetAddress(&m_MC_imot1);

    if ((b = m_tree_MC->FindBranch("MC_vx")))
        b->SetAddress(&m_MC_vx);

    if ((b = m_tree_MC->FindBranch("MC_vy")))
        b->SetAddress(&m_MC_vy);

    if ((b = m_tree_MC->FindBranch("MC_vz")))
        b->SetAddress(&m_MC_vz);
    
    
    // Set up redundant branches, which should be dropped eventually
    if ((b = m_tree_MC->FindBranch("MC_index")))
        b->SetAddress(&m_MC_index);

    if ((b = m_tree_MC->FindBranch("MC_e")))
        b->SetAddress(&m_MC_E);

    if ((b = m_tree_MC->FindBranch("MC_px")))
        b->SetAddress(&m_MC_px);

    if ((b = m_tree_MC->FindBranch("MC_py")))
        b->SetAddress(&m_MC_py);

    if ((b = m_tree_MC->FindBranch("MC_pz")))
        b->SetAddress(&m_MC_pz);

    if ((b = m_tree_MC->FindBranch("MC_eta")))
        b->SetAddress(&m_MC_eta);

    if ((b = m_tree_MC->FindBranch("MC_phi")))
        b->SetAddress(&m_MC_phi);

    if (m_doJpsi)
    {
        if ((b = m_tree_MC->FindBranch("MC_JPsiFromTop")))
            b->SetAddress(&m_MC_JPsiFromTop);
      
        if ((b = m_tree_MC->FindBranch("MC_JPsiFromAntiTop")))
            b->SetAddress(&m_MC_JPsiFromAntiTop);
      
        if ((b = m_tree_MC->FindBranch("MC_LeptonFromTop")))
            b->SetAddress(&m_MC_LeptonFromTop);
      
        if ((b = m_tree_MC->FindBranch("MC_LeptonFromAntiTop")))
            b->SetAddress(&m_MC_LeptonFromAntiTop);
    }
    
    
    // Since some branches are optional and might not be available, fill the corresponding arrays
    //with defaults
    for (unsigned i = 0; i < unsigned(m_n_MCs); ++i)
    {
        m_MC_status[i] = 0;
        m_MC_isLastCopy[i] = false;
        m_MC_imot1[i] = -1;
        m_MC_vx[i] = 0.f;
        m_MC_vy[i] = 0.f;
        m_MC_vz[i] = 0.f;
        
        // Redundant arrays
        m_MC_index[i] = i;  // see a comment on this array in the writeInfo method
        m_MC_E[i] = 0.f;
        m_MC_px[i] = 0.f;
        m_MC_py[i] = 0.f;
        m_MC_pz[i] = 0.f;
        m_MC_eta[i] = 0.f;
        m_MC_phi[i] = 0.f;
        m_MC_LeptonFromTop[i] = false;
        m_MC_LeptonFromAntiTop[i] = false;
        m_MC_JPsiFromTop[i] = false;
        m_MC_JPsiFromAntiTop[i] = false;
    }
}


MCExtractor::~MCExtractor()
{
    delete m_MC_lorentzvector;
}


void MCExtractor::doConsumes(edm::ConsumesCollector&& collector)
{
    SuperBaseExtractor::doConsumes(std::forward<edm::ConsumesCollector>(collector));
    
    m_genParticleToken = collector.consumes<edm::View<reco::GenParticle>>(m_genParticleTag);
}


void MCExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor*)
{
    // Read the collection of pruned generator-level particles
    edm::Handle<edm::View<reco::GenParticle>> genParticlesHandle;
    event.getByToken(m_genParticleToken, genParticlesHandle);
    edm::View<reco::GenParticle> const &genParticles = *genParticlesHandle;
    
    
    // (Absolute) PDG IDs of some particles to avoid using magic numbers in the code
    static const int pdgIdE = 11, pdgIdMu = 13, pdgIdPhoton = 22, pdgIdTop = 6, pdgIdJpsi = 443;
    

    // Loop over particles and choose a small subset to store. Save them in an intermediate
    //collection instead of the output tree in order to setup properly parental relations
    vector<reco::GenParticle const *> selectedParticles;
    
    for (auto const &p: genParticles)
    {
        bool selected = false;
        int const absPdgId = abs(p.pdgId());
        
        
        // Particles from the final state of the matrix element (the definition is generator-
        //dependent)
        if (p.fromHardProcessBeforeFSR())
            selected = true;
        
        // Prompt electrons, muons, or photons
        if (p.isPromptFinalState() and
         (absPdgId == pdgIdE or absPdgId == pdgIdMu or absPdgId == pdgIdPhoton))
            selected = true;
        
        if (m_doJpsi)
        {
            // J/psi decaying into a pair of muons
            if (absPdgId == pdgIdJpsi and p.numberOfDaughters() == 2 and
             abs(p.daughter(0)->pdgId()) == pdgIdMu and abs(p.daughter(0)->pdgId()) == pdgIdMu)
                selected = true;
            
            // Muon from the decay of a J/psi
            if (absPdgId == pdgIdMu and p.numberOfMothers() > 0 and
             p.mother(0)->pdgId() == pdgIdJpsi)
                selected = true;
        }
        
        
        if (selected)
            selectedParticles.emplace_back(&p);
    }
    
    
    // Save basic information about selected particles in the output tree
    m_n_MCs = std::min<unsigned>(selectedParticles.size(), m_MCs_MAX);
    m_MC_lorentzvector->Clear();
    
    for (unsigned i = 0; i < unsigned(m_n_MCs); ++i)
    {
        auto const &p = *selectedParticles.at(i);
        
        m_MC_type[i] = p.pdgId();
        m_MC_status[i] = p.status();
        m_MC_isLastCopy[i] = p.isLastCopy();
        new((*m_MC_lorentzvector)[i]) TLorentzVector(p.px(), p.py(), p.pz(), p.energy());
        m_MC_vx[i] = p.vx();
        m_MC_vy[i] = p.vy();
        m_MC_vz[i] = p.vz();
        
        // Duplicates. Provided for backward compatibility only. Should be removed in future
        m_MC_E[i] = p.energy();
        m_MC_px[i] = p.px();
        m_MC_py[i] = p.py();
        m_MC_pz[i] = p.pz();
        m_MC_eta[i] = p.eta();
        m_MC_phi[i] = p.phi();
    }
    
    
    // For each selected particle, find its closest ancestor among other selected particles and
    //store its position in the vector selectedParticles. Since the relative ordering of selected
    //particles is same as in the original collection in AOD, ancestors of a particle are located
    //closer to the beginning of the vector than the particle. This also means that the first
    //particle in the vector has no stored ancestors
    m_MC_imot1[0] = -1;
    
    for (unsigned iPart = unsigned(m_n_MCs) - 1; iPart > 0; --iPart)
    {
        reco::Candidate const *m = selectedParticles.at(iPart);
        bool motherFound = false;
        
        
        // Check the full inheritance line of the current particle
        while (not motherFound)
        {
            // Move one step back in history if possible
            if (m->numberOfMothers() == 0)
                break;
            
            m = m->mother(0);
            
            
            // Check if the current ancestor is stored in the vector
            for (unsigned iMotherCand = iPart; iMotherCand-- > 0; )
            {
                if (m == selectedParticles.at(iMotherCand))
                {
                    m_MC_imot1[iPart] = iMotherCand;
                    motherFound = true;
                    break;
                }
            }
        }
        
        
        // None of the selected particles is an ancestor of the current particle
        if (not motherFound)
            m_MC_imot1[iPart] = -1;
    }
    
    
    // In the past, mothers were identified by their indices in the original collection in AOD.
    //Consequently, these indices had to be stored for all selected particles. Now mother indices
    //correspond to actual arrays saved in the output trees, but for the sake of backward
    //compatibility identifying indices are kept. They should be removed in future
    for (unsigned i = 0; i < unsigned(m_n_MCs); ++i)
        m_MC_index[i] = i;
    
    
    // Fill several arrays specific for an analysis with J/psi. The same information can be obtained
    //from stored ancestors. It is only saved for backward compatibility and should be removed in
    //future
    if (m_doJpsi)
    {
        // Set all flags to false by default
        for (unsigned i = 0; i < unsigned(m_n_MCs); ++i)
        {
            m_MC_LeptonFromTop[i] = false;
            m_MC_LeptonFromAntiTop[i] = false;
            m_MC_JPsiFromTop[i] = false;
            m_MC_JPsiFromAntiTop[i] = false;
        }
        
        
        // Loop over selected particles
        for (unsigned i = 0; i < unsigned(m_n_MCs); ++i)
        {
            bool isLepton = false, isJpsi = false;
            int const absPdgId = abs(m_MC_type[i]);
            
            
            // Only consider leptons and J/psi
            if (absPdgId == pdgIdE or absPdgId == pdgIdMu)
                isLepton = true;
            else if (absPdgId == pdgIdJpsi)
                isJpsi = true;
            else
                continue;
            
            
            // Check if the current particle originates from a top quark or antiquark
            int ancestorIndex = i;
            
            while (true)
            {
                // Move one ancestor back if possible
                ancestorIndex = m_MC_imot1[ancestorIndex];
                
                if (ancestorIndex == -1)
                    break;
                
                
                // Check if the ancestor is a top quark and update the flags accordingly
                int ancestorPdgId = selectedParticles.at(ancestorIndex)->pdgId();
                
                if (abs(ancestorPdgId) == pdgIdTop)
                {
                    if (isLepton)
                    {
                        if (ancestorPdgId > 0)
                            m_MC_LeptonFromTop[i] = true;
                        else
                            m_MC_LeptonFromAntiTop[i] = true;
                    }
                    else if (isJpsi)
                    {
                        if (ancestorPdgId > 0)
                            m_MC_JPsiFromTop[i] = true;
                        else
                            m_MC_JPsiFromAntiTop[i] = true;
                    }
                    
                    break;
                }
            }
        }
    }
    
        
    // Debug print out
    #if 0
    cout << "\n\n=== New event ===\n";
    
    for (unsigned i = 0; i < unsigned(m_n_MCs); ++i)
    {
        cout << " #" << i << ": PDG ID: " << m_MC_type[i] << ", status: " << m_MC_status[i] << '\n';
        
        TLorentzVector const &p4 = *dynamic_cast<TLorentzVector *>(m_MC_lorentzvector->At(i));
        cout << "  pt: " << p4.Pt();
        
        if (p4.Pt() > 0.)
            cout << ", eta: " << p4.Eta() << ", phi: " << p4.Phi();
        
        cout << "\n  first mother: " << m_MC_imot1[i] << '\n';
        cout << "  last copy? " << m_MC_isLastCopy[i] << '\n';
        cout << "  J/psi flags: " << m_MC_LeptonFromTop[i] << ", " << m_MC_LeptonFromAntiTop[i] <<
         ", " << m_MC_JPsiFromTop[i] << ", " << m_MC_JPsiFromAntiTop[i] << endl;
    }
    #endif
    
    
    // Finally, write the arrays in the output tree
    fillTree();
}


void MCExtractor::getInfo(int ievt) 
{
    m_tree_MC->GetEntry(ievt); 
}


void MCExtractor::fillTree()
{
  m_tree_MC->Fill(); 
}


int MCExtractor::getSize() const
{
    return m_n_MCs;
}


int MCExtractor::getType(int index) const
{
    return m_MC_type[index];
}


int MCExtractor::getStatus(int index) const
{
    return m_MC_status[index];
}


bool MCExtractor::getIsLastCopy(int index) const
{
    return m_MC_isLastCopy[index];
}


TLorentzVector const &MCExtractor::getP4(int index) const
{
    return *dynamic_cast<TLorentzVector*>((*m_MC_lorentzvector)[index]);
}


int MCExtractor::getMom1Index(int index) const
{
    return m_MC_imot1[index];
}


int MCExtractor::getPatIndex(int index) const
{
    return m_MC_index[index];
}


float MCExtractor::getPx(int index) const
{
    return getP4(index).Px();
}


float MCExtractor::getPy(int index) const
{
    return getP4(index).Py();
}


float MCExtractor::getPz(int index) const
{
    return getP4(index).Pz();
}


float MCExtractor::getE(int index) const
{
    return getP4(index).E();
}
