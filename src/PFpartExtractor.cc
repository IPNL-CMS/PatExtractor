#include "../interface/PFpartExtractor.h"


PFpartExtractor::PFpartExtractor(const std::string& name, const edm::InputTag& tag, bool doTree)
  : BaseExtractor(name)
{
  m_tag = tag;

  // Set everything to 0
  m_OK = false;
  m_pf_lorentzvector = new TClonesArray("TLorentzVector");
  reset();

  // Tree definition

  if (doTree)
  {
    m_OK = true;
    m_tree_pfpart         = new TTree(name.c_str(), "PF particles info");     
    m_tree_pfpart->Branch("n_pf",      &m_pf_size,   "n_pf/I");  
    m_tree_pfpart->Branch("pf_4vector","TClonesArray",&m_pf_lorentzvector, 100, 0);
    m_tree_pfpart->Branch("pf_vx",     &m_pf_vx,     "pf_vx[n_pf]/F");  
    m_tree_pfpart->Branch("pf_vy",     &m_pf_vy,     "pf_vy[n_pf]/F");  
    m_tree_pfpart->Branch("pf_vz",     &m_pf_vz,     "pf_vz[n_pf]/F");  
    m_tree_pfpart->Branch("pf_charge", &m_pf_charge, "pf_charge[n_pf]/I");
    m_tree_pfpart->Branch("pf_pdgid",  &m_pf_pdgid,  "pf_charge[n_pf]/I");

  }
}

PFpartExtractor::PFpartExtractor(const std::string& name, TFile *a_file)
  :BaseExtractor(name)
{
  std::cout << "PFpartExtractor objet is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_pfpart = dynamic_cast<TTree*>(a_file->Get(m_name.c_str()));

  if (!m_tree_pfpart)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_pf_lorentzvector = new TClonesArray("TLorentzVector");

  m_tree_pfpart->SetBranchAddress("n_pfs",     &m_pf_size);
  m_tree_pfpart->SetBranchAddress("pf_4vector",&m_pf_lorentzvector);
  m_tree_pfpart->SetBranchAddress("pf_vx",     &m_pf_vx);
  m_tree_pfpart->SetBranchAddress("pf_vy",     &m_pf_vy);
  m_tree_pfpart->SetBranchAddress("pf_vz",     &m_pf_vz);
  m_tree_pfpart->SetBranchAddress("pf_charge", &m_pf_charge);
  m_tree_pfpart->SetBranchAddress("pf_pdgid",  &m_pf_pdgid);

}

PFpartExtractor::~PFpartExtractor()
{}



//
// Method filling the main particle tree
//

void PFpartExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const reco::PFCandidate& part, int index) 
{
  if ( abs(part.pdgId()) != 13 ) return;

  if (m_pf_size>=m_pfpart_MAX) return;

  new((*m_pf_lorentzvector)[m_pf_size]) TLorentzVector(part.px(),part.py(),part.pz(),part.energy());
  if ( !part.trackRef() ) {
    m_pf_vx[index]              = 0;
    m_pf_vy[index]              = 0;
    m_pf_vz[index]              = 0;
  } else {
    m_pf_vx[index]              = part.trackRef()->vx();
    m_pf_vy[index]              = part.trackRef()->vy();
    m_pf_vz[index]              = part.trackRef()->vz();
  }
  m_pf_charge[m_pf_size]          = part.charge();
  m_pf_pdgid[m_pf_size]           = part.pdgId();

  ++m_pf_size;

}


// Method initializing everything (to do for each event)

void PFpartExtractor::reset()
{
  m_size = 0;
  m_pf_size = 0;

  for (int i=0;i<m_pfpart_MAX;++i) 
  {
    m_pf_vx[i]         = 0.;
    m_pf_vy[i]         = 0.;
    m_pf_vz[i]         = 0.;
    m_pf_charge[i]     = 0;
    m_pf_pdgid[i]      = 0.;
  }
  
  m_pf_lorentzvector->Clear();

}


void PFpartExtractor::fillTree()
{
  m_tree_pfpart->Fill(); 
}
 
//
// Method getting the info from an input file
//

void PFpartExtractor::getInfo(int ievt) 
{
  m_tree_pfpart->GetEntry(ievt); 
}
