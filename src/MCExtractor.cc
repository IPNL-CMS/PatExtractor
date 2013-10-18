#include "../interface/MCExtractor.h"


MCExtractor::MCExtractor(const std::string& name, bool doTree, bool doJpsi)
{
  _doJpsi = doJpsi; 
  // Set everything to 0

  m_OK = false;
  m_MC_lorentzvector = new TClonesArray("TLorentzVector");
  reset();

  // Tree definition

  if (doTree)
  {
    m_OK = true;

    m_tree_MC = new TTree(name.c_str(), "PAT MC info");  
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
    if (_doJpsi) {
      m_tree_MC->Branch("MC_JPsiFromTop",  &m_MC_JPsiFromTop,  "m_MC_JPsiFromTop[n_MCs]/B");  
      m_tree_MC->Branch("MC_JPsiFromAntiTop",  &m_MC_JPsiFromAntiTop,  "m_MC_JPsiFromAntiTop[n_MCs]/B");  
      m_tree_MC->Branch("MC_LeptonFromTop",  &m_MC_LeptonFromTop,  "m_MC_LeptonFromTop[n_MCs]/B");  
      m_tree_MC->Branch("MC_LeptonFromAntiTop",  &m_MC_LeptonFromAntiTop,  "m_MC_LeptonFromAntiTop[n_MCs]/B");
    }
  }
}

MCExtractor::MCExtractor(const std::string& name, TFile *a_file, bool doJpsi)
{
  _doJpsi = doJpsi; 
  std::cout << "MCExtractor object is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_MC = dynamic_cast<TTree*>(a_file->Get(name.c_str()));

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
  if (_doJpsi) {
    if (m_tree_MC->FindBranch("MC_JPsiFromTop")) 
      m_tree_MC->SetBranchAddress("MC_JPsiFromTop", &m_MC_JPsiFromTop);
    if (m_tree_MC->FindBranch("MC_JPsiFromAntiTop")) 
      m_tree_MC->SetBranchAddress("MC_JPsiFromAntiTop", &m_MC_JPsiFromAntiTop);
    if (m_tree_MC->FindBranch("MC_LeptonFromTop")) 
      m_tree_MC->SetBranchAddress("MC_LeptonFromTop", &m_MC_LeptonFromTop);
    if (m_tree_MC->FindBranch("MC_LeptonFromAntiTop")) 
      m_tree_MC->SetBranchAddress("MC_LeptonFromAntiTop", &m_MC_LeptonFromAntiTop);
  }
}


MCExtractor::~MCExtractor()
{}



//
// Method filling the main particle tree
//


void MCExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor)
{
  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByLabel("genParticles", genParticles);

  MCExtractor::reset();
  m_n_MCs = static_cast<int>(genParticles->size());

  int   id = 0;
  int   id_r = 0;
  float px_r = 0.;
  float py_r = 0.;
  float pz_r = 0.;
  int st    = 0;
  int n_mot = 0;
  //int n_dau = 0;

  int ipart = 0;

  const reco::GenParticle* mothertmp;
  const reco::GenParticle* motherleptontmp;

  for(int i=0; i < m_n_MCs; ++i) 
  {
    const reco::Candidate & p = (*genParticles)[i];

    id    = p.pdgId();
    st    = p.status(); 
    n_mot = p.numberOfMothers(); 
    //n_dau = p.numberOfDaughters();

    int iMo1 = -1;
    int iMo2 = -1;

    if (st == 3 || ( _doJpsi && id == 443 && p.numberOfDaughters() == 2 && fabs(p.daughter(0)->pdgId()) ==  13 && fabs(p.daughter(1)->pdgId()) ==  13))
    {

      // MC@NLO use different status code
      // Status 3 are only assigned to proton, with Px = Py = Pt = 0.
      // Remove then until proper support of MC@NLO
      if (p.px() == 0 && p.py() == 0)
        continue;

      if (n_mot > 0)
      {
        id_r = (p.mother(0))->pdgId();
        px_r = (p.mother(0))->px();
        py_r = (p.mother(0))->py();
        pz_r = (p.mother(0))->pz();

        for(int j = 0; j < m_n_MCs; ++j) 
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

          for(int j=0; j < m_n_MCs; ++j) 
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
      new((*m_MC_lorentzvector)[ipart]) TLorentzVector(p.px(),p.py(),p.pz(),p.energy());

      if (_doJpsi && id==443) {
        mothertmp = &(*genParticles)[i];
        for (int ifrag=0; ifrag<10; ifrag++) {
          if (mothertmp->mother() == 0) break;
          mothertmp = (reco::GenParticle*) mothertmp->mother();
          if (fabs(mothertmp->pdgId())==92 || fabs(mothertmp->pdgId())==91) break;
        }
        for (int imb=0; imb<100; imb++) {
          if (mothertmp->mother(imb) != 0 && fabs(mothertmp->mother(imb)->pdgId())==5) {
            mothertmp = (reco::GenParticle*) mothertmp->mother(imb);
            break;
          } 
        }
        if (fabs(mothertmp->pdgId()) == 5) {
          for (int ib=0; ib<10; ib++) {
            if (mothertmp->mother() !=0 && (fabs(mothertmp->mother()->pdgId())==5 || fabs(mothertmp->mother()->pdgId())==21)) mothertmp = (reco::GenParticle*) mothertmp->mother();
            else break;
          }
          if(mothertmp->mother() !=0) mothertmp = (reco::GenParticle*) mothertmp->mother();
          if (mothertmp->pdgId()==6)  m_MC_JPsiFromTop[ipart] = true;
          if (mothertmp->pdgId()==-6) m_MC_JPsiFromAntiTop[ipart] = true;
        }	
        //std::cout << "JPsiFromTop = " << m_MC_JPsiFromTop[ipart] << "JPsiFromAntiTop = " << m_MC_JPsiFromAntiTop[ipart] << std::endl;
      }	
      if (_doJpsi && (fabs(id)==11 || fabs(id) ==13)) {
        motherleptontmp = &(*genParticles)[i];
        for (int i=0; i<10; i++) {
          if (motherleptontmp->mother() == 0) break;
          motherleptontmp = (reco::GenParticle*) motherleptontmp->mother();
          if (motherleptontmp->pdgId()==6) {
            m_MC_LeptonFromTop[ipart] = true;
            break;
          }
          if (motherleptontmp->pdgId()==-6) {
            m_MC_LeptonFromAntiTop[ipart] = true;
            break;
          }
        }
        //std::cout << "LeptonFromTop = " << m_MC_LeptonFromTop[ipart] << "LeptonFromAntiTop = " << m_MC_LeptonFromAntiTop[ipart] << std::endl;
        if (m_MC_LeptonFromTop[ipart]==false && m_MC_LeptonFromAntiTop[ipart]==false) std::cout << "cas suspect" << std::endl;
      }


      if (n_mot==0) m_MC_generation[ipart] = 0;

      ++ipart;
    }
  }  


  m_n_MCs = ipart;

  for(int i=1; i<6; ++i) 
    MCExtractor::constructGeneration(i, m_n_MCs);

  fillTree();
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
    if (_doJpsi) {
      m_MC_JPsiFromTop[i] = false;
      m_MC_JPsiFromAntiTop[i] = false;
      m_MC_LeptonFromTop[i] = false;
      m_MC_LeptonFromAntiTop[i] = false;
    }
  }
  m_MC_lorentzvector->Clear();
}


void MCExtractor::fillTree()
{
  m_tree_MC->Fill(); 
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

int MCExtractor::getSize() const {
  return m_n_MCs;
}
