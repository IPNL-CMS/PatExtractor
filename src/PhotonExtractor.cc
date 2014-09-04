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
    m_tree_photon->Branch("photon_hasPixelSeed",  &m_pho_hasPixelSeed, "photon_hasPixelSeed[n_photons]/O");
    m_tree_photon->Branch("photon_hadTowOverEm",  &m_pho_hadTowOverEm,   "photon_hadTowOverEm[n_photons]/F");
    m_tree_photon->Branch("photon_sigmaIetaIeta",  &m_pho_sigmaIetaIeta,   "photon_sigmaIetaIeta[n_photons]/F");
    m_tree_photon->Branch("photon_hasMatchedPromptElectron",  &m_pho_hasMatchedPromptElectron,   "photon_hasMatchedPromptElectron[n_photons]/O");
    m_tree_photon->Branch("photon_chargedHadronsIsolation", &m_pho_chargedHadronsIsolation, "photon_chargedHadronsIsolation[n_photons]/F");
    m_tree_photon->Branch("photon_neutralHadronsIsolation",  &m_pho_neutralHadronsIsolation,   "photon_neutralHadronsIsolation[n_photons]/F");
    m_tree_photon->Branch("photon_photonIsolation",  &m_pho_photonIsolation,   "photon_photonIsolation[n_photons]/F");
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
  m_tree_photon->SetBranchAddress("photon_hasPixelSeed",  &m_pho_hasPixelSeed);
  m_tree_photon->SetBranchAddress("photon_hadTowOverEm",  &m_pho_hadTowOverEm);
  m_tree_photon->SetBranchAddress("photon_sigmaIetaIeta",  &m_pho_sigmaIetaIeta);
  m_tree_photon->SetBranchAddress("photon_hasMatchedPromptElectron",  &m_pho_hasMatchedPromptElectron);
  m_tree_photon->SetBranchAddress("photon_chargedHadronsIsolation",  &m_pho_chargedHadronsIsolation);
  m_tree_photon->SetBranchAddress("photon_neutralHadronsIsolation",  &m_pho_neutralHadronsIsolation);
  m_tree_photon->SetBranchAddress("photon_photonIsolation",  &m_pho_photonIsolation);
  m_tree_photon->SetBranchAddress("photon_mcParticleIndex",&m_pho_MCIndex);  
}


PhotonExtractor::~PhotonExtractor()
{}


void PhotonExtractor::doConsumes(edm::ConsumesCollector&& collector) {
  BaseExtractor::doConsumes(std::forward<edm::ConsumesCollector>(collector));

  m_rhoToken = collector.consumes<double>(edm::InputTag("kt6PFJets", "rho", "RECO"));
  m_matchedPromptElectronToken = collector.consumes<edm::ValueMap<bool>>(edm::InputTag("photonPFIsolation", "hasMatchedPromptElectron", "PAT"));
  m_chargedHadronsIsolationToken = collector.consumes<edm::ValueMap<double>>(edm::InputTag("photonPFIsolation", "chargedHadronsIsolation", "PAT"));
  m_neutralHadronsIsolationToken = collector.consumes<edm::ValueMap<double>>(edm::InputTag("photonPFIsolation", "neutralHadronsIsolation", "PAT"));
  m_photonsIsolationToken = collector.consumes<edm::ValueMap<double>>(edm::InputTag("photonPFIsolation", "photonsIsolation", "PAT"));
}

enum class IsolationType {
  CHARGED_HADRONS,
  NEUTRAL_HADRONS,
  PHOTONS
};

float getEffectiveArea(float eta, IsolationType type) {
  eta = fabs(eta);
  switch (type) {
    case IsolationType::CHARGED_HADRONS:
      if (eta < 1.0)
        return 0.012;
      else if (eta < 1.479)
        return 0.010;
      else if (eta < 2.0)
        return 0.014;
      else if (eta < 2.2)
        return 0.012;
      else if (eta < 2.3)
        return 0.016;
      else if (eta < 2.4)
        return 0.020;
      else
        return 0.012;
      break;

    case IsolationType::NEUTRAL_HADRONS:
      if (eta < 1.0)
        return 0.030;
      else if (eta < 1.479)
        return 0.057;
      else if (eta < 2.0)
        return 0.039;
      else if (eta < 2.2)
        return 0.015;
      else if (eta < 2.3)
        return 0.024;
      else if (eta < 2.4)
        return 0.039;
      else
        return 0.072;
      break;

    case IsolationType::PHOTONS:
      if (eta < 1.0)
        return 0.148;
      else if (eta < 1.479)
        return 0.130;
      else if (eta < 2.0)
        return 0.112;
      else if (eta < 2.2)
        return 0.216;
      else if (eta < 2.3)
        return 0.262;
      else if (eta < 2.4)
        return 0.260;
      else
        return 0.266;
      break;
  }

  return -1;
}

double getCorrectedPFIsolation(double isolation, double rho, float eta, IsolationType type) {
  float effectiveArea = getEffectiveArea(eta, type);

  return std::max(isolation - rho * effectiveArea, 0.);
}




//
// Method filling the main particle tree
//

void PhotonExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* m_MC) 
{
  edm::Handle<pat::PhotonCollection>  photonHandle;
  event.getByToken(m_token, photonHandle);
  pat::PhotonCollection p_photons = *photonHandle;

  reset();
  m_size = 0;

  //std::cout<<"write info photon before loop"<<std::endl;
  for (unsigned int i = 0; i < p_photons.size(); ++i)
  {
    //std::cout<<"photon #"<<i<<std::endl;
    pat::PhotonRef photonRef(photonHandle, i);    
    PhotonExtractor::writeInfo(event, iSetup, p_photons.at(i), m_size, photonRef); 

    if (m_MC)
      doMCMatch(p_photons.at(i), event, m_MC, m_size);

    m_size++;
  }

  fillTree();
}

void PhotonExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Photon& part, int index, const pat::PhotonRef& photonRef) 
{
  if (index>=m_photons_MAX) return;
  //std::cout<<"write photon #"<<index<<std::endl;

  new((*m_pho_lorentzvector)[index]) TLorentzVector(part.px(),part.py(),part.pz(),part.energy());
  m_pho_vx[index]   = part.vx();
  m_pho_vy[index]   = part.vy();
  m_pho_vz[index]   = part.vz();
  m_pho_hasPixelSeed[index] = part.hasPixelSeed();
  m_pho_hadTowOverEm[index]  = photonRef->hadTowOverEm();
  m_pho_sigmaIetaIeta[index] = photonRef->sigmaIetaIeta();
  
  edm::Handle<double> rhos;
  event.getByToken(m_rhoToken, rhos);
  double rho = *rhos;

  // Isolations are produced at PAT level by the P\u1e27otonPFIsolation producer
  edm::Handle<edm::ValueMap<bool>> hasMatchedPromptElectronHandle;
  event.getByToken(m_matchedPromptElectronToken, hasMatchedPromptElectronHandle);
  m_pho_hasMatchedPromptElectron[index] = (*hasMatchedPromptElectronHandle)[photonRef];

  // Now, isolations
  edm::Handle<edm::ValueMap<double>> chargedHadronsIsolationHandle;
  event.getByToken(m_chargedHadronsIsolationToken, chargedHadronsIsolationHandle);

  edm::Handle<edm::ValueMap<double>> neutralHadronsIsolationHandle;
  event.getByToken(m_neutralHadronsIsolationToken, neutralHadronsIsolationHandle);

  edm::Handle<edm::ValueMap<double>> photonIsolationHandle;
  event.getByToken(m_photonsIsolationToken, photonIsolationHandle);

  m_pho_chargedHadronsIsolation[index] = getCorrectedPFIsolation((*chargedHadronsIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::CHARGED_HADRONS);
  m_pho_neutralHadronsIsolation[index] = getCorrectedPFIsolation((*neutralHadronsIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS);
  m_pho_photonIsolation[index]         = getCorrectedPFIsolation((*photonIsolationHandle)[photonRef], rho, photonRef->eta(), IsolationType::PHOTONS);

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
    m_pho_hasPixelSeed[i] = false;
    m_pho_hadTowOverEm[i] = -1;
    m_pho_sigmaIetaIeta[i] = -1;
    m_pho_hasMatchedPromptElectron[i] = false;
    m_pho_chargedHadronsIsolation[i] = -1;
    m_pho_neutralHadronsIsolation[i] = -1;
    m_pho_photonIsolation[i] = -1;
    m_pho_MCIndex[i] = -1;
  }
  
  if(m_pho_lorentzvector)
    m_pho_lorentzvector->Clear();
}


void PhotonExtractor::fillTree()
{
  m_tree_photon->Fill(); 
}
