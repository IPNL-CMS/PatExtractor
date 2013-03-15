#include "../interface/PhotonExtractor.h"


PhotonExtractor::PhotonExtractor(const std::string& name, const edm::InputTag& tag, bool doTree)
  : BaseExtractor(name)
{
  m_tag = tag;
  m_deltaR_cut = 0.2; // Maximum acceptable distance for MC matching
  
  // Set everything to 0

  m_pho_lorentzvector = new TClonesArray("TLorentzVector");
  reset();


  // Tree definition

  if (doTree)
  {
    m_tree_photon   = new TTree(name.c_str(), "PAT photon info");  
    m_tree_photon->Branch("n_photons",  &m_size, "n_photons/i");  
    m_tree_photon->Branch("photon_4vector","TClonesArray",&m_pho_lorentzvector, 1000, 0);    
    m_tree_photon->Branch("photon_vx",  &m_pho_vx,   "photon_vx[n_photons]/F");  
    m_tree_photon->Branch("photon_vy",  &m_pho_vy,   "photon_vy[n_photons]/F");  
    m_tree_photon->Branch("photon_vz",  &m_pho_vz,   "photon_vz[n_photons]/F");
    m_tree_photon->Branch("photon_mcParticleIndex",&m_pho_MCIndex,"photon_mcParticleIndex[n_photons]/I");  
  }
}


PhotonExtractor::PhotonExtractor(const std::string& name, TFile *a_file)
  : BaseExtractor(name)
{
  std::cout << "PhotonExtractor objet is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_photon = dynamic_cast<TTree*>(a_file->Get(name.c_str()));

  if (!m_tree_photon)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  } 

  m_OK = true;

  m_pho_lorentzvector = new TClonesArray("TLorentzVector");


  // Branches definition

  m_tree_photon->SetBranchAddress("n_photons",  &m_size);
  m_tree_photon->SetBranchAddress("photon_4vector",&m_pho_lorentzvector);
  m_tree_photon->SetBranchAddress("photon_vx",  &m_pho_vx);
  m_tree_photon->SetBranchAddress("photon_vy",  &m_pho_vy);
  m_tree_photon->SetBranchAddress("photon_vz",  &m_pho_vz);
  m_tree_photon->SetBranchAddress("photon_mcParticleIndex",&m_pho_MCIndex);  
}


PhotonExtractor::~PhotonExtractor()
{}



//
// Method filling the main particle tree
//

void PhotonExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Photon& part, int index) 
{
  if (index>=m_photons_MAX) return;

  new((*m_pho_lorentzvector)[index]) TLorentzVector(part.energy(),part.px(),part.py(),part.pz());
  m_pho_vx[index]   = part.vx();
  m_pho_vy[index]   = part.vy();
  m_pho_vz[index]   = part.vz();
}


//
// Method getting the info from an input file
//

void PhotonExtractor::getInfo(int ievt) 
{
  m_tree_photon->GetEntry(ievt); 
}

// Method initializing everything (to do for each event)

void PhotonExtractor::reset()
{
  m_size = 0;

  for (int i=0;i<m_photons_MAX;++i) 
  {
    m_pho_vx[i] = 0.;
    m_pho_vy[i] = 0.;
    m_pho_vz[i] = 0.;
    m_pho_MCIndex[i] = -1;
  }
  m_pho_lorentzvector->Clear();
}


void PhotonExtractor::fillTree()
{
  m_tree_photon->Fill(); 
}
