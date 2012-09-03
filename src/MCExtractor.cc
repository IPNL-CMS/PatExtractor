#include "../interface/MCExtractor.h"


MCExtractor::MCExtractor(bool doTree)
{
  // Set everything to 0

  m_OK = false;
  m_MC_lorentzvector = new TClonesArray("TLorentzVector");
  MCExtractor::reset();


  // Tree definition

  if (doTree)
  {
    m_OK = true;

    m_tree_MC = new TTree("MC","PAT MC info");  
    m_tree_MC->Branch("MC_4vector","TClonesArray",&m_MC_lorentzvector, 1000, 0);
    m_tree_MC->Branch("n_MCs",  &m_n_MCs,"n_MCs/I");  
    m_tree_MC->Branch("MC_index",   &m_MC_index,    "MC_index[n_MCs]/I");  
    m_tree_MC->Branch("MC_type",    &m_MC_type,     "MC_type[n_MCs]/I");  
    m_tree_MC->Branch("MC_mot1",    &m_MC_imot1,    "MC_mot1[n_MCs]/I");  
    m_tree_MC->Branch("MC_mot2",    &m_MC_imot2,    "MC_mot2[n_MCs]/I");  
    m_tree_MC->Branch("MC_generation",   &m_MC_generation,    "MC_generation[n_MCs]/I");  
    m_tree_MC->Branch("MC_e",   &m_MC_E,    "MC_e[n_MCs]/F");  
    m_tree_MC->Branch("MC_px",  &m_MC_px,   "MC_px[n_MCs]/F");  
    m_tree_MC->Branch("MC_py",  &m_MC_py,   "MC_py[n_MCs]/F");  
    m_tree_MC->Branch("MC_pz",  &m_MC_pz,   "MC_pz[n_MCs]/F");  
    m_tree_MC->Branch("MC_vx",  &m_MC_vx,   "MC_vx[n_MCs]/F");  
    m_tree_MC->Branch("MC_vy",  &m_MC_vy,   "MC_vy[n_MCs]/F");  
    m_tree_MC->Branch("MC_vz",  &m_MC_vz,   "MC_vz[n_MCs]/F");
    m_tree_MC->Branch("MC_eta", &m_MC_eta,  "MC_eta[n_MCs]/F");  
    m_tree_MC->Branch("MC_phi", &m_MC_phi,  "MC_phi[n_MCs]/F");  
  }
}

MCExtractor::MCExtractor(TFile *a_file)
{
  std::cout << "MCExtractor object is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_MC = dynamic_cast<TTree*>(a_file->Get("MC"));

  if (!m_tree_MC)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }
  
  m_OK = true;

  m_MC_lorentzvector = new TClonesArray("TLorentzVector");

  if (m_tree_MC->FindBranch("n_MCs")) 
    m_tree_MC->SetBranchAddress("n_MCs",  &m_n_MCs);
  if (m_tree_MC->FindBranch("MC_4vector")) 
    m_tree_MC->SetBranchAddress("MC_4vector",&m_MC_lorentzvector);
  if (m_tree_MC->FindBranch("MC_index")) 
    m_tree_MC->SetBranchAddress("MC_index",   &m_MC_index);
  if (m_tree_MC->FindBranch("MC_type")) 
    m_tree_MC->SetBranchAddress("MC_type",    &m_MC_type);
  if (m_tree_MC->FindBranch("MC_mot1")) 
    m_tree_MC->SetBranchAddress("MC_mot1",    &m_MC_imot1);
  if (m_tree_MC->FindBranch("MC_mot2")) 
    m_tree_MC->SetBranchAddress("MC_mot2",    &m_MC_imot2);
  if (m_tree_MC->FindBranch("MC_generation")) 
    m_tree_MC->SetBranchAddress("MC_generation",   &m_MC_generation);
  if (m_tree_MC->FindBranch("MC_e")) 
    m_tree_MC->SetBranchAddress("MC_e",   &m_MC_E);
  if (m_tree_MC->FindBranch("MC_px")) 
    m_tree_MC->SetBranchAddress("MC_px",  &m_MC_px);
  if (m_tree_MC->FindBranch("MC_py")) 
    m_tree_MC->SetBranchAddress("MC_py",  &m_MC_py);
  if (m_tree_MC->FindBranch("MC_pz")) 
    m_tree_MC->SetBranchAddress("MC_pz",  &m_MC_pz);
  if (m_tree_MC->FindBranch("MC_vx")) 
    m_tree_MC->SetBranchAddress("MC_vx",  &m_MC_vx);
  if (m_tree_MC->FindBranch("MC_vy")) 
    m_tree_MC->SetBranchAddress("MC_vy",  &m_MC_vy);
  if (m_tree_MC->FindBranch("MC_vz")) 
    m_tree_MC->SetBranchAddress("MC_vz",  &m_MC_vz);
  if (m_tree_MC->FindBranch("MC_eta")) 
    m_tree_MC->SetBranchAddress("MC_eta", &m_MC_eta);
  if (m_tree_MC->FindBranch("MC_phi")) 
    m_tree_MC->SetBranchAddress("MC_phi", &m_MC_phi);
}


MCExtractor::~MCExtractor()
{}



//
// Method filling the main particle tree
//


void MCExtractor::writeInfo(const edm::Event *event) 
{
  edm::Handle<reco::GenParticleCollection> genParticles;
  event->getByLabel("genParticles", genParticles);
  
  MCExtractor::reset();
  MCExtractor::fillSize(static_cast<int>(genParticles->size()));

  int   id = 0;
  int   id_r = 0;
  float px_r = 0.;
  float py_r = 0.;
  float pz_r = 0.;
  int st    = 0;
  int n_mot = 0;
  int n_dau = 0;

  int ipart = 0;

  for(int i=0; i<MCExtractor::getSize(); ++i) 
  {
    const reco::Candidate & p = (*genParticles)[i];
    
    id    = p.pdgId();
    st    = p.status(); 
    n_mot = p.numberOfMothers(); 
    n_dau = p.numberOfDaughters();
    
    int iMo1 = -1;
    int iMo2 = -1;
    
    if (st==3)
    {
      if (n_mot>0)
      {
	id_r = (p.mother(0))->pdgId();
	px_r = (p.mother(0))->px();
	py_r = (p.mother(0))->py();
	pz_r = (p.mother(0))->pz();
	
	for(int j=0; j<m_n_MCs; ++j) 
        {
	  const reco::Candidate &p2 = (*genParticles)[j];

	  if (p2.pdgId() != id_r) continue;
	  if (fabs(p2.px()-px_r)>0.0001) continue;
	  if (fabs(p2.py()-py_r)>0.0001) continue;
	  if (fabs(p2.pz()-pz_r)>0.0001) continue;
	  
	  iMo1=j;
	  
	  break;
	}

	if (n_mot>1)
        {
	  id_r = (p.mother(1))->pdgId();
	  px_r = (p.mother(1))->px();
	  py_r = (p.mother(1))->py();
	  pz_r = (p.mother(1))->pz();
	  
	  for(int j=0; j<MCExtractor::getSize(); ++j) 
          {
	    const reco::Candidate &p2 = (*genParticles)[j];
	    
	    if (p2.pdgId() != id_r) continue;
	    if (fabs(p2.px()-px_r)>0.0001) continue;
	    if (fabs(p2.py()-py_r)>0.0001) continue;
	    if (fabs(p2.pz()-pz_r)>0.0001) continue;
	    
	    iMo2=j;
	    
	    break;
	  }
	}
      }
      
      m_MC_imot1[ipart]      = iMo1;
      m_MC_imot2[ipart]      = iMo2;
      m_MC_index[ipart]      = i;
      m_MC_status[ipart]     = st;
      m_MC_type[ipart]       = id;
      m_MC_E[ipart]          = p.energy();
      m_MC_px[ipart]         = p.px();
      m_MC_py[ipart]         = p.py();
      m_MC_pz[ipart]         = p.pz();
      m_MC_vx[ipart]         = p.vx();
      m_MC_vy[ipart]         = p.vy();
      m_MC_vz[ipart]         = p.vz();
      m_MC_eta[ipart]        = p.eta();
      m_MC_phi[ipart]        = p.phi();
      new((*m_MC_lorentzvector)[i]) TLorentzVector(p.px(),p.py(),p.pz(),p.energy());

      if (n_mot==0) m_MC_generation[ipart] = 0;
      
      ++ipart;
    }
  }  
  

  MCExtractor::fillSize(ipart);
  
  for(int i=1; i<6; ++i) 
  {    
    MCExtractor::constructGeneration(i,MCExtractor::getSize());
  }
  
  MCExtractor::fillTree();
 
}

//
// Method getting the info from an input file
//

void MCExtractor::getInfo(int ievt) 
{
  m_tree_MC->GetEntry(ievt); 
}

// Method initializing everything (to do for each event)

void MCExtractor::reset()
{
  m_n_MCs = 0;
  
  for (int i=0;i<m_MCs_MAX;++i) 
  {
    m_MC_index[i] = 0;
    m_MC_status[i] = 0;
    m_MC_type[i] = 0;
    m_MC_imot1[i] = 0;
    m_MC_imot2[i] = 0;
    m_MC_generation[i] = -1;
    m_MC_E[i] = 0.;
    m_MC_px[i] = 0.;
    m_MC_py[i] = 0.;
    m_MC_pz[i] = 0.;
    m_MC_vx[i] = 0.;
    m_MC_vy[i] = 0.;
    m_MC_vz[i] = 0.;
    m_MC_eta[i] = 0.;
    m_MC_phi[i] = 0.; 
  }
  m_MC_lorentzvector->Clear();
}


void MCExtractor::fillTree()
{
  m_tree_MC->Fill(); 
}
 
void MCExtractor::fillSize(int size)
{
  m_n_MCs=size;
}

int  MCExtractor::getSize()
{
  return m_n_MCs;
}


void MCExtractor::constructGeneration(int gene, int npart)
{
  for(int i=0; i<npart; ++i) 
  {
    if (m_MC_generation[i]==gene-1)
    {
      int index = m_MC_index[i];

      for(int j=0; j<npart; ++j) 
      {
        if (m_MC_imot1[j]==index) m_MC_generation[j]=gene;
        if (m_MC_imot2[j]==index) m_MC_generation[j]=gene;
      }
    }
  }
}
