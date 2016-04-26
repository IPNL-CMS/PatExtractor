#include "../interface/MCExtractor.h"


MCExtractor::MCExtractor(const std::string& name, bool doTree, bool doJpsi, bool doD0)
{
  _doJpsi = doJpsi; 
  _doD0 = doD0;
  // Set everything to 0

  m_OK = false;
  m_MC_lorentzvector = new TClonesArray("TLorentzVector");
  if (_doJpsi){
    m_MC_JPsi_lorentzvector = new TClonesArray("TLorentzVector");
    m_MC_Bhad_lorentzvector = new TClonesArray("TLorentzVector");
    m_MC_Bquark_lorentzvector = new TClonesArray("TLorentzVector");
  }
  if (_doD0) {
    m_MC_D0_lorentzvector = new TClonesArray("TLorentzVector");
    m_MC_D0_daughter0_lorentzvector = new TClonesArray("TLorentzVector");
    m_MC_D0_daughter1_lorentzvector = new TClonesArray("TLorentzVector");
    m_MC_D0_Bhad_lorentzvector = new TClonesArray("TLorentzVector");
  }
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
      m_tree_MC->Branch("MC_JPsi_4vector","TClonesArray",&m_MC_JPsi_lorentzvector, 1000, 0);
      m_tree_MC->Branch("MC_Bhad_4vector","TClonesArray",&m_MC_Bhad_lorentzvector, 1000, 0);
      m_tree_MC->Branch("MC_Bhad_id",  &m_MC_Bhad_id,   "MC_Bhad_id[n_MCs]/I");  
      m_tree_MC->Branch("MC_BhadWithNuDaughter",  &m_MC_BhadWithNuDaughter, "MC_BhadWithNuDaughter[n_MCs]/B");  
      m_tree_MC->Branch("MC_BhadWithoutNuDaughter",  &m_MC_BhadWithoutNuDaughter, "MC_BhadWithoutNuDaughter[n_MCs]/B");  
      m_tree_MC->Branch("MC_Bquark_4vector","TClonesArray",&m_MC_Bquark_lorentzvector, 1000, 0);
    }
    if (_doD0) {
      m_tree_MC->Branch("MC_D0_4vector","TClonesArray",&m_MC_D0_lorentzvector, 1000, 0);
      m_tree_MC->Branch("MC_D0_daughter0_4vector","TClonesArray",&m_MC_D0_daughter0_lorentzvector, 1000, 0);
      m_tree_MC->Branch("MC_D0_daughter0_id",  &m_MC_D0_daughter0_id,   "MC_daughter0_id[n_MCs]/I");  
      m_tree_MC->Branch("MC_D0_daughter1_4vector","TClonesArray",&m_MC_D0_daughter1_lorentzvector, 1000, 0);
      m_tree_MC->Branch("MC_D0_daughter1_id",  &m_MC_D0_daughter1_id,   "MC_D0_daughter1_id[n_MCs]/I");  
      m_tree_MC->Branch("MC_D0_Bhad_4vector","TClonesArray",&m_MC_D0_Bhad_lorentzvector, 1000, 0);
      m_tree_MC->Branch("MC_D0_Bhad_id",  &m_MC_D0_Bhad_id,   "MC_Bhad_id[n_MCs]/I");  
    }
  }

  idBHadrons.push_back(411); idBHadrons.push_back(421); idBHadrons.push_back(10411); idBHadrons.push_back(10421); 
  idBHadrons.push_back(413); idBHadrons.push_back(423); idBHadrons.push_back(10413); idBHadrons.push_back(10423); 
  idBHadrons.push_back(20413); idBHadrons.push_back(20423); idBHadrons.push_back(415); idBHadrons.push_back(425); 
  idBHadrons.push_back(431); idBHadrons.push_back(10431); idBHadrons.push_back(433); idBHadrons.push_back(10433); 
  idBHadrons.push_back(20433); idBHadrons.push_back(435); 
  //bottom mesons
  idBHadrons.push_back(511); idBHadrons.push_back(521); idBHadrons.push_back(10511); idBHadrons.push_back(10521); 
  idBHadrons.push_back(513); idBHadrons.push_back(523); idBHadrons.push_back(10513); idBHadrons.push_back(10523); 
  idBHadrons.push_back(20513); idBHadrons.push_back(20523); idBHadrons.push_back(515); idBHadrons.push_back(525); 
  idBHadrons.push_back(531); idBHadrons.push_back(10531); idBHadrons.push_back(533); idBHadrons.push_back(10533); 
  idBHadrons.push_back(20533); idBHadrons.push_back(535); idBHadrons.push_back(541); idBHadrons.push_back(10541); 
  idBHadrons.push_back(543); idBHadrons.push_back(10543); idBHadrons.push_back(20543); idBHadrons.push_back(545); 
  // ccbar mesons
  idBHadrons.push_back(441); idBHadrons.push_back(10441); idBHadrons.push_back(100441); idBHadrons.push_back(443); 
  idBHadrons.push_back(10443); idBHadrons.push_back(20443); idBHadrons.push_back(100443); idBHadrons.push_back(30443); 
  idBHadrons.push_back(9000443); idBHadrons.push_back(9010443); idBHadrons.push_back(9020443); idBHadrons.push_back(445); 
  idBHadrons.push_back(100445); 
  //bbar mesons:w
  idBHadrons.push_back(551); idBHadrons.push_back(100551); idBHadrons.push_back(110551); idBHadrons.push_back(200551); 
  idBHadrons.push_back(210551); idBHadrons.push_back(553); idBHadrons.push_back(10553); idBHadrons.push_back(20553); 
  idBHadrons.push_back(30553); idBHadrons.push_back(100553); idBHadrons.push_back(110553); idBHadrons.push_back(120553); 
  idBHadrons.push_back(130553); idBHadrons.push_back(200553); idBHadrons.push_back(210553); idBHadrons.push_back(220553); 
  idBHadrons.push_back(300553); idBHadrons.push_back(9000553); idBHadrons.push_back(9010553); idBHadrons.push_back(555); 
  idBHadrons.push_back(10555); idBHadrons.push_back(20555); idBHadrons.push_back(100555); idBHadrons.push_back(110555); 
  idBHadrons.push_back(120555); idBHadrons.push_back(200555); idBHadrons.push_back(557); idBHadrons.push_back(100557); 
  // charmed baryons
  idBHadrons.push_back(4122); idBHadrons.push_back(4222); idBHadrons.push_back(4212); idBHadrons.push_back(4112); 
  idBHadrons.push_back(4224); idBHadrons.push_back(4214); idBHadrons.push_back(4114); idBHadrons.push_back(4232); 
  idBHadrons.push_back(4132); idBHadrons.push_back(4322); idBHadrons.push_back(4312); idBHadrons.push_back(4324); 
  idBHadrons.push_back(4314); idBHadrons.push_back(4332); idBHadrons.push_back(4334); idBHadrons.push_back(4412); 
  idBHadrons.push_back(4422); idBHadrons.push_back(4414); idBHadrons.push_back(4424); idBHadrons.push_back(4432); 
  idBHadrons.push_back(4434); idBHadrons.push_back(4444);  
  // bottom baryons
  idBHadrons.push_back(5122); idBHadrons.push_back(5112); idBHadrons.push_back(5212); idBHadrons.push_back(5222); 
  idBHadrons.push_back(5114); idBHadrons.push_back(5214); idBHadrons.push_back(5224); idBHadrons.push_back(5132); 
  idBHadrons.push_back(5232); idBHadrons.push_back(5312); idBHadrons.push_back(5322); idBHadrons.push_back(5314); 
  idBHadrons.push_back(5324); idBHadrons.push_back(5332); idBHadrons.push_back(5334); idBHadrons.push_back(5142); 
  idBHadrons.push_back(5242); idBHadrons.push_back(5412); idBHadrons.push_back(5422); idBHadrons.push_back(5414); 
  idBHadrons.push_back(5424); idBHadrons.push_back(5342); idBHadrons.push_back(5432); idBHadrons.push_back(5434); 
  idBHadrons.push_back(5442); idBHadrons.push_back(5444); idBHadrons.push_back(5512); idBHadrons.push_back(5522); 
  idBHadrons.push_back(5514); idBHadrons.push_back(5524); idBHadrons.push_back(5532); idBHadrons.push_back(5534); 
  idBHadrons.push_back(5542); idBHadrons.push_back(5544); idBHadrons.push_back(5554);  

}

MCExtractor::MCExtractor(const std::string& name, TFile *a_file, bool doJpsi, bool doD0)
{
  _doJpsi = doJpsi; 
  _doD0 = doD0; 
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
  if(_doJpsi) {
    m_MC_JPsi_lorentzvector = new TClonesArray("TLorentzVector");
    m_MC_Bhad_lorentzvector = new TClonesArray("TLorentzVector");
    m_MC_Bquark_lorentzvector = new TClonesArray("TLorentzVector");
  }
  if(_doD0) {
    m_MC_D0_lorentzvector = new TClonesArray("TLorentzVector");
    m_MC_D0_daughter0_lorentzvector = new TClonesArray("TLorentzVector");
    m_MC_D0_daughter1_lorentzvector = new TClonesArray("TLorentzVector");
    m_MC_D0_Bhad_lorentzvector = new TClonesArray("TLorentzVector");
  }

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
    if (m_tree_MC->FindBranch("MC_JPsi_4vector")) 
      m_tree_MC->SetBranchAddress("MC_JPsi_4vector",&m_MC_JPsi_lorentzvector);
    if (m_tree_MC->FindBranch("MC_Bhad_4vector")) 
      m_tree_MC->SetBranchAddress("MC_Bhad_4vector",&m_MC_Bhad_lorentzvector);
    if (m_tree_MC->FindBranch("MC_Bhad_id")) 
      m_tree_MC->SetBranchAddress("MC_Bhad_id",  &m_MC_Bhad_id);
    if (m_tree_MC->FindBranch("MC_BhadWithNuDaughter"))
      m_tree_MC->SetBranchAddress("MC_BhadWithNuDaughter",  &m_MC_BhadWithNuDaughter);  
    if (m_tree_MC->FindBranch("MC_BhadWithoutNuDaughter"))
      m_tree_MC->SetBranchAddress("MC_BhadWithoutNuDaughter",  &m_MC_BhadWithoutNuDaughter);  
    if (m_tree_MC->FindBranch("MC_Bquark_4vector")) 
      m_tree_MC->SetBranchAddress("MC_Bquark_4vector",&m_MC_Bquark_lorentzvector);
  }
  if (_doD0) {
    if (m_tree_MC->FindBranch("MC_D0_4vector")) 
      m_tree_MC->SetBranchAddress("MC_D0_4vector",&m_MC_D0_lorentzvector);
    if (m_tree_MC->FindBranch("MC_D0_daughter0_4vector")) 
      m_tree_MC->SetBranchAddress("MC_D0_daughter0_4vector",&m_MC_D0_daughter0_lorentzvector);
    if (m_tree_MC->FindBranch("MC_D0_daughter0_id")) 
      m_tree_MC->SetBranchAddress("MC_D0_daughter0_id",  &m_MC_D0_daughter0_id);
    if (m_tree_MC->FindBranch("MC_D0_daughter1_4vector")) 
      m_tree_MC->SetBranchAddress("MC_D0_daughter1_4vector",&m_MC_D0_daughter1_lorentzvector);
    if (m_tree_MC->FindBranch("MC_D0_daughter1_id")) 
      m_tree_MC->SetBranchAddress("MC_D0_daughter1_id",  &m_MC_D0_daughter1_id);
    if (m_tree_MC->FindBranch("MC_D0_Bhad_4vector")) 
      m_tree_MC->SetBranchAddress("MC_D0_Bhad_4vector",&m_MC_D0_Bhad_lorentzvector);
    if (m_tree_MC->FindBranch("MC_D0_Bhad_id")) 
      m_tree_MC->SetBranchAddress("MC_D0_Bhad_id",  &m_MC_D0_Bhad_id);
  }

  idBHadrons.push_back(411); idBHadrons.push_back(421); idBHadrons.push_back(10411); idBHadrons.push_back(10421); 
  idBHadrons.push_back(413); idBHadrons.push_back(423); idBHadrons.push_back(10413); idBHadrons.push_back(10423); 
  idBHadrons.push_back(20413); idBHadrons.push_back(20423); idBHadrons.push_back(415); idBHadrons.push_back(425); 
  idBHadrons.push_back(431); idBHadrons.push_back(10431); idBHadrons.push_back(433); idBHadrons.push_back(10433); 
  idBHadrons.push_back(20433); idBHadrons.push_back(435); 
  //bottom mesons
  idBHadrons.push_back(511); idBHadrons.push_back(521); idBHadrons.push_back(10511); idBHadrons.push_back(10521); 
  idBHadrons.push_back(513); idBHadrons.push_back(523); idBHadrons.push_back(10513); idBHadrons.push_back(10523); 
  idBHadrons.push_back(20513); idBHadrons.push_back(20523); idBHadrons.push_back(515); idBHadrons.push_back(525); 
  idBHadrons.push_back(531); idBHadrons.push_back(10531); idBHadrons.push_back(533); idBHadrons.push_back(10533); 
  idBHadrons.push_back(20533); idBHadrons.push_back(535); idBHadrons.push_back(541); idBHadrons.push_back(10541); 
  idBHadrons.push_back(543); idBHadrons.push_back(10543); idBHadrons.push_back(20543); idBHadrons.push_back(545); 
  // ccbar mesons
  idBHadrons.push_back(441); idBHadrons.push_back(10441); idBHadrons.push_back(100441); idBHadrons.push_back(443); 
  idBHadrons.push_back(10443); idBHadrons.push_back(20443); idBHadrons.push_back(100443); idBHadrons.push_back(30443); 
  idBHadrons.push_back(9000443); idBHadrons.push_back(9010443); idBHadrons.push_back(9020443); idBHadrons.push_back(445); 
  idBHadrons.push_back(100445); 
  //bbar mesons:w
  idBHadrons.push_back(551); idBHadrons.push_back(100551); idBHadrons.push_back(110551); idBHadrons.push_back(200551); 
  idBHadrons.push_back(210551); idBHadrons.push_back(553); idBHadrons.push_back(10553); idBHadrons.push_back(20553); 
  idBHadrons.push_back(30553); idBHadrons.push_back(100553); idBHadrons.push_back(110553); idBHadrons.push_back(120553); 
  idBHadrons.push_back(130553); idBHadrons.push_back(200553); idBHadrons.push_back(210553); idBHadrons.push_back(220553); 
  idBHadrons.push_back(300553); idBHadrons.push_back(9000553); idBHadrons.push_back(9010553); idBHadrons.push_back(555); 
  idBHadrons.push_back(10555); idBHadrons.push_back(20555); idBHadrons.push_back(100555); idBHadrons.push_back(110555); 
  idBHadrons.push_back(120555); idBHadrons.push_back(200555); idBHadrons.push_back(557); idBHadrons.push_back(100557); 
  // charmed baryons
  idBHadrons.push_back(4122); idBHadrons.push_back(4222); idBHadrons.push_back(4212); idBHadrons.push_back(4112); 
  idBHadrons.push_back(4224); idBHadrons.push_back(4214); idBHadrons.push_back(4114); idBHadrons.push_back(4232); 
  idBHadrons.push_back(4132); idBHadrons.push_back(4322); idBHadrons.push_back(4312); idBHadrons.push_back(4324); 
  idBHadrons.push_back(4314); idBHadrons.push_back(4332); idBHadrons.push_back(4334); idBHadrons.push_back(4412); 
  idBHadrons.push_back(4422); idBHadrons.push_back(4414); idBHadrons.push_back(4424); idBHadrons.push_back(4432); 
  idBHadrons.push_back(4434); idBHadrons.push_back(4444);  
  // bottom baryons
  idBHadrons.push_back(5122); idBHadrons.push_back(5112); idBHadrons.push_back(5212); idBHadrons.push_back(5222); 
  idBHadrons.push_back(5114); idBHadrons.push_back(5214); idBHadrons.push_back(5224); idBHadrons.push_back(5132); 
  idBHadrons.push_back(5232); idBHadrons.push_back(5312); idBHadrons.push_back(5322); idBHadrons.push_back(5314); 
  idBHadrons.push_back(5324); idBHadrons.push_back(5332); idBHadrons.push_back(5334); idBHadrons.push_back(5142); 
  idBHadrons.push_back(5242); idBHadrons.push_back(5412); idBHadrons.push_back(5422); idBHadrons.push_back(5414); 
  idBHadrons.push_back(5424); idBHadrons.push_back(5342); idBHadrons.push_back(5432); idBHadrons.push_back(5434); 
  idBHadrons.push_back(5442); idBHadrons.push_back(5444); idBHadrons.push_back(5512); idBHadrons.push_back(5522); 
  idBHadrons.push_back(5514); idBHadrons.push_back(5524); idBHadrons.push_back(5532); idBHadrons.push_back(5534); 
  idBHadrons.push_back(5542); idBHadrons.push_back(5544); idBHadrons.push_back(5554);  

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

    bool isBhad = false;
    for (unsigned int iBHad = 0; iBHad < idBHadrons.size(); iBHad++) {
      if (abs(id) != idBHadrons[iBHad]) continue;
      bool willRadiate = false;
      for (unsigned int imb = 0; imb < p.numberOfDaughters(); imb++) {
        if (p.daughter(imb) != 0 && abs(p.daughter(imb)->pdgId()) == 22) {
          willRadiate = true;
          break;
        }
      }
      if (!willRadiate) isBhad = true;
    }

    if (st == 3 || (st >= 21 && st <= 29) || (_doJpsi && id == 443 && p.numberOfDaughters() == 2 && abs(p.daughter(0)->pdgId()) ==  13 && abs(p.daughter(1)->pdgId()) ==  13) || (_doD0 && id == 421 && p.numberOfDaughters() == 2 && ((abs(p.daughter(0)->pdgId()) ==  321 && abs(p.daughter(1)->pdgId()) ==  211) || (abs(p.daughter(0)->pdgId()) ==  211 && abs(p.daughter(1)->pdgId()) ==  321))) || isBhad)
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

      if (_doD0 && id == 421) {
        mothertmp = &(*genParticles)[i];
        new((*m_MC_D0_lorentzvector)[ipart]) TLorentzVector(mothertmp->px(),mothertmp->py(),mothertmp->pz(),mothertmp->energy());
        new((*m_MC_D0_daughter0_lorentzvector)[ipart]) TLorentzVector(mothertmp->daughter(0)->px(),mothertmp->daughter(0)->py(),mothertmp->daughter(0)->pz(),mothertmp->daughter(0)->energy());
        m_MC_D0_daughter0_id[ipart] = mothertmp->daughter(0)->pdgId();
        new((*m_MC_D0_daughter1_lorentzvector)[ipart]) TLorentzVector(mothertmp->daughter(1)->px(),mothertmp->daughter(1)->py(),mothertmp->daughter(1)->pz(),mothertmp->daughter(1)->energy());
        m_MC_D0_daughter1_id[ipart] = mothertmp->daughter(1)->pdgId();
        for (int ifrag=0; ifrag<10; ifrag++) {
          if (mothertmp->mother() == 0) break;
          mothertmp = (reco::GenParticle*) mothertmp->mother();
          if (ifrag==0) {
            new((*m_MC_Bhad_lorentzvector)[ipart]) TLorentzVector(mothertmp->px(),mothertmp->py(),mothertmp->pz(),mothertmp->energy());
            m_MC_Bhad_id[ipart] = mothertmp->pdgId();
            for (unsigned int idau = 0; idau < mothertmp->numberOfDaughters(); idau++) {
              if (isNeutrinoPdgId(mothertmp->daughter(idau)->pdgId())) m_MC_BhadWithNuDaughter[ipart] = true;
            }
            if (!m_MC_BhadWithNuDaughter) m_MC_BhadWithoutNuDaughter[ipart] = true;
          }
          if (mothertmp->mother() != 0 && (abs(mothertmp->mother()->pdgId())==92 || abs(mothertmp->mother()->pdgId())==91)) {
            new((*m_MC_D0_Bhad_lorentzvector)[ipart]) TLorentzVector(mothertmp->px(),mothertmp->py(),mothertmp->pz(),mothertmp->energy());
            m_MC_D0_Bhad_id[ipart] = mothertmp->pdgId();
            break;
          }
        }
      }
      if (_doJpsi && id == 443) {
        mothertmp = &(*genParticles)[i];
        new((*m_MC_JPsi_lorentzvector)[ipart]) TLorentzVector(mothertmp->px(),mothertmp->py(),mothertmp->pz(),mothertmp->energy());
        for (int ifrag=0; ifrag<10; ifrag++) {
          if (mothertmp->mother() == 0) break;
          mothertmp = (reco::GenParticle*) mothertmp->mother();
          if (ifrag==0) {
            new((*m_MC_Bhad_lorentzvector)[ipart]) TLorentzVector(mothertmp->px(),mothertmp->py(),mothertmp->pz(),mothertmp->energy());
            m_MC_Bhad_id[ipart] = mothertmp->pdgId();
            for (unsigned int idau = 0; idau < mothertmp->numberOfDaughters(); idau++) {
              if (isNeutrinoPdgId(mothertmp->daughter(idau)->pdgId())) m_MC_BhadWithNuDaughter[ipart] = true;
            }
            if (!m_MC_BhadWithNuDaughter) m_MC_BhadWithoutNuDaughter[ipart] = true;
          }
          if (abs(mothertmp->pdgId())==92 || abs(mothertmp->pdgId())==91) break;
        }
        for (int imb=0; imb<100; imb++) {
          if (mothertmp->mother(imb) != 0 && abs(mothertmp->mother(imb)->pdgId())==5) {
            mothertmp = (reco::GenParticle*) mothertmp->mother(imb);
            new((*m_MC_Bquark_lorentzvector)[ipart]) TLorentzVector(mothertmp->px(),mothertmp->py(),mothertmp->pz(),mothertmp->energy());
            break;
          } 
        }
        if (abs(mothertmp->pdgId()) == 5) {
          for (int ib=0; ib<10; ib++) {
            if (mothertmp->mother() !=0 && (abs(mothertmp->mother()->pdgId())==5 || abs(mothertmp->mother()->pdgId())==21)) mothertmp = (reco::GenParticle*) mothertmp->mother();
            else break;
          }
          if(mothertmp->mother() !=0) mothertmp = (reco::GenParticle*) mothertmp->mother();
          if (mothertmp->pdgId()==6)  m_MC_JPsiFromTop[ipart] = true;
          if (mothertmp->pdgId()==-6) m_MC_JPsiFromAntiTop[ipart] = true;
        }	
        //std::cout << "JPsiFromTop = " << m_MC_JPsiFromTop[ipart] << "JPsiFromAntiTop = " << m_MC_JPsiFromAntiTop[ipart] << std::endl;
      }	
      if (_doJpsi && (abs(id) == 11 || abs(id) == 13)) {
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
      m_MC_Bhad_id[i] = 0;
      m_MC_BhadWithNuDaughter[i] = false;
      m_MC_BhadWithoutNuDaughter[i] = false;
    }
    if (_doD0) {
      m_MC_D0_daughter0_id[i] = 0;
      m_MC_D0_daughter1_id[i] = 0;
      m_MC_D0_Bhad_id[i] = 0;
    }
  }
  m_MC_lorentzvector->Clear();
  if (_doJpsi) {
    m_MC_JPsi_lorentzvector->Clear();
    m_MC_Bhad_lorentzvector->Clear();
    m_MC_Bquark_lorentzvector->Clear();
  }
  if (_doD0) {
    m_MC_D0_lorentzvector->Clear();
    m_MC_D0_daughter0_lorentzvector->Clear();
    m_MC_D0_daughter1_lorentzvector->Clear();
    m_MC_D0_Bhad_lorentzvector->Clear();
  }
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
