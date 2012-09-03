#include "../interface/PhotonExtractor.h"


PhotonExtractor::PhotonExtractor(bool doTree, edm::InputTag tag)
{
  m_tag = tag;
  m_deltaR_cut = 0.2; // Maximum acceptable distance for MC matching
  
  // Set everything to 0

  m_pho_lorentzvector = new TClonesArray("TLorentzVector");
  PhotonExtractor::reset();


  // Tree definition

  if (doTree)
  {
    m_tree_photon   = new TTree("photon","PAT photon info");  
    m_tree_photon->Branch("n_photons",  &m_n_photons,"n_photons/I");  
    m_tree_photon->Branch("photon_4vector","TClonesArray",&m_pho_lorentzvector, 1000, 0);    
    m_tree_photon->Branch("photon_vx",  &m_pho_vx,   "photon_vx[n_photons]/F");  
    m_tree_photon->Branch("photon_vy",  &m_pho_vy,   "photon_vy[n_photons]/F");  
    m_tree_photon->Branch("photon_vz",  &m_pho_vz,   "photon_vz[n_photons]/F");
    m_tree_photon->Branch("photon_mcParticleIndex",&m_pho_MCIndex,"photon_mcParticleIndex[n_photons]/I");  
  }
}


PhotonExtractor::PhotonExtractor(TFile *a_file)
{
  std::cout << "PhotonExtractor objet is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_photon = dynamic_cast<TTree*>(a_file->Get("photon"));

  if (!m_tree_photon)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  } 

  m_OK = true;

  m_pho_lorentzvector = new TClonesArray("TLorentzVector");


  // Branches definition

  m_tree_photon->SetBranchAddress("n_photons",  &m_n_photons);
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

void PhotonExtractor::writeInfo(const edm::Event *event,MCExtractor* m_MC, bool doMC) 
{
  edm::Handle< edm::View<pat::Photon> >  photonHandle;
  event->getByLabel(m_tag, photonHandle);
  edm::View<pat::Photon> p_photons = *photonHandle;

  PhotonExtractor::reset();
  PhotonExtractor::fillSize(static_cast<int>(p_photons.size()));

  if (PhotonExtractor::getSize())
  {
    for(int i=0; i<PhotonExtractor::getSize(); ++i)
    {
      PhotonExtractor::writeInfo(&p_photons.at(i),i);
      
      if (doMC)
      {
	int idx_min = PhotonExtractor::getMatch(&p_photons.at(i),m_MC); 
	m_pho_MCIndex[i] = idx_min;
      }
    }
  }

  PhotonExtractor::fillTree();
}


void PhotonExtractor::writeInfo(const pat::Photon *part, int index) 
{
  if (index>=m_photons_MAX) return;

  new((*m_pho_lorentzvector)[index]) TLorentzVector(part->energy(),part->px(),part->py(),part->pz());
  m_pho_vx[index]   = part->vx();
  m_pho_vy[index]   = part->vy();
  m_pho_vz[index]   = part->vz();
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
  m_n_photons = 0;

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
 
void PhotonExtractor::fillSize(int size)
{
  m_n_photons=size;
}

int  PhotonExtractor::getSize()
{
  return m_n_photons;
}

int PhotonExtractor::getMatch(const pat::Photon *part, MCExtractor* m_MC)
{
  float deltaR_min = 1e6;
  int idx_min    = -1;
      
  for(int mcPart_i=0; mcPart_i<m_MC->getSize(); ++mcPart_i) 
  {
    if (m_MC->getStatus(mcPart_i)!=3) continue;
    if (fabs(m_MC->getType(mcPart_i))!=22) continue;
    

    TLorentzVector TL_genPart(m_MC->getPx(mcPart_i),m_MC->getPy(mcPart_i),m_MC->getPz(mcPart_i),m_MC->getE(mcPart_i));
    TLorentzVector TL_photon(part->px(),part->py(),part->pz(),part->energy());
    	   
    if(TL_genPart.Pt())
    {
      float deltaR = TL_genPart.DeltaR(TL_photon);
      //float deltaP = fabs(TL_genPart.Pt()-TL_photon.Pt());

      if(deltaR<deltaR_min)
      {
	deltaR_min = deltaR;
	idx_min = mcPart_i;
      }
    }
  }
  
  if (deltaR_min>m_deltaR_cut)
    return -2;

  return idx_min;

}
