#include "../interface/PhotonExtractor.h"


PhotonExtractor::PhotonExtractor(const std::string& name, const edm::ParameterSet& settings)
  : BaseExtractor(name, settings),
  m_photonTag(settings.getParameter<edm::InputTag>("input")),
  m_rhoTag(settings.getParameter<edm::InputTag>("rho")),
  m_phoLooseIdMapTag(settings.getParameter<edm::InputTag>("phoLooseIdMap")),
  m_phoMediumIdMapTag(settings.getParameter<edm::InputTag>("phoMediumIdMap")),
  m_phoTightIdMapTag(settings.getParameter<edm::InputTag>("phoTightIdMap"))
{
  m_deltaR_cut = 0.2; // Maximum acceptable distance for MC matching
  
  // Set everything to 0

  m_pho_lorentzvector = new TClonesArray("TLorentzVector");
  reset();


  // Tree definition
  m_tree_photon   = new TTree(name.c_str(), "PAT photon info");  
  m_tree_photon->Branch("n_photons",  &m_size, "n_photons/i");  
  m_tree_photon->Branch("photon_4vector","TClonesArray",&m_pho_lorentzvector, 1000, 0);    
  m_tree_photon->Branch("photon_vx",  &m_pho_vx,   "photon_vx[n_photons]/F");  
  m_tree_photon->Branch("photon_vy",  &m_pho_vy,   "photon_vy[n_photons]/F");  
  m_tree_photon->Branch("photon_vz",  &m_pho_vz,   "photon_vz[n_photons]/F");
  m_tree_photon->Branch("photon_passLooseId",  &m_pho_passLooseId, "photon_passLooseId[n_photons]/O");
  m_tree_photon->Branch("photon_passMediumId",  &m_pho_passMediumId, "photon_passMediumId[n_photons]/O");
  m_tree_photon->Branch("photon_passTightId",  &m_pho_passTightId, "photon_passTightId[n_photons]/O");
  m_tree_photon->Branch("photon_chargedHadronsIsolation", &m_pho_chargedHadronsIsolation, "photon_chargedHadronsIsolation[n_photons]/F");
  m_tree_photon->Branch("photon_neutralHadronsIsolation",  &m_pho_neutralHadronsIsolation,   "photon_neutralHadronsIsolation[n_photons]/F");
  m_tree_photon->Branch("photon_photonIsolation",  &m_pho_photonIsolation,   "photon_photonIsolation[n_photons]/F");
  m_tree_photon->Branch("photon_mcParticleIndex",&m_pho_MCIndex,"photon_mcParticleIndex[n_photons]/I");  
  
  
  // Mark that the extractor has been constructed properly
  setHealthy(true);
}


PhotonExtractor::PhotonExtractor(const std::string& name, const edm::ParameterSet& settings, TFile *a_file)
  : BaseExtractor(name, settings)
{
  std::cout << "PhotonExtractor objet is retrieved" << std::endl;

  // Tree definition
  setHealthy(false);
  m_tree_photon = dynamic_cast<TTree*>(a_file->Get(name.c_str()));

  if (!m_tree_photon)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  } 

  m_pho_lorentzvector = new TClonesArray("TLorentzVector");


  // Branches definition

  m_tree_photon->SetBranchAddress("n_photons",  &m_size);
  m_tree_photon->SetBranchAddress("photon_4vector",&m_pho_lorentzvector);
  m_tree_photon->SetBranchAddress("photon_vx",  &m_pho_vx);
  m_tree_photon->SetBranchAddress("photon_vy",  &m_pho_vy);
  m_tree_photon->SetBranchAddress("photon_vz",  &m_pho_vz);
  m_tree_photon->SetBranchAddress("photon_passLooseId",  &m_pho_passLooseId);
  m_tree_photon->SetBranchAddress("photon_passMediumId",  &m_pho_passMediumId);
  m_tree_photon->SetBranchAddress("photon_passTightId",  &m_pho_passTightId);
  m_tree_photon->SetBranchAddress("photon_chargedHadronsIsolation",  &m_pho_chargedHadronsIsolation);
  m_tree_photon->SetBranchAddress("photon_neutralHadronsIsolation",  &m_pho_neutralHadronsIsolation);
  m_tree_photon->SetBranchAddress("photon_photonIsolation",  &m_pho_photonIsolation);
  m_tree_photon->SetBranchAddress("photon_mcParticleIndex",&m_pho_MCIndex);
  
  
  // Mark that the extractor has been constructed properly
  setHealthy(true);
}


PhotonExtractor::~PhotonExtractor()
{}


void PhotonExtractor::doConsumes(edm::ConsumesCollector&& collector) {
  BaseExtractor::doConsumes(std::forward<edm::ConsumesCollector>(collector));

  m_token = collector.consumes<pat::PhotonCollection>(m_photonTag);
  m_rhoToken = collector.consumes<double>(m_rhoTag);
  m_phoLooseIdMapToken = collector.consumes<edm::ValueMap<bool> >(m_phoLooseIdMapTag);
  m_phoMediumIdMapToken = collector.consumes<edm::ValueMap<bool> >(m_phoMediumIdMapTag);
  m_phoTightIdMapToken = collector.consumes<edm::ValueMap<bool> >(m_phoTightIdMapTag);
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
        return 0.0234;
      else if (eta < 1.479)
        return 0.0189;
      else if (eta < 2.0)
        return 0.0171;
      else if (eta < 2.2)
        return 0.0129;
      else if (eta < 2.3)
        return 0.0110;
      else if (eta < 2.4)
        return 0.0074;
      else
        return 0.0035;
      break;

    case IsolationType::NEUTRAL_HADRONS:
      if (eta < 1.0)
        return 0.0053;
      else if (eta < 1.479)
        return 0.0103;
      else if (eta < 2.0)
        return 0.0057;
      else if (eta < 2.2)
        return 0.0070;
      else if (eta < 2.3)
        return 0.0152;
      else if (eta < 2.4)
        return 0.0232;
      else
        return 0.1709;
      break;

    case IsolationType::PHOTONS:
      if (eta < 1.0)
        return 0.078;
      else if (eta < 1.479)
        return 0.0629;
      else if (eta < 2.0)
        return 0.0264;
      else if (eta < 2.2)
        return 0.0462;
      else if (eta < 2.3)
        return 0.0740;
      else if (eta < 2.4)
        return 0.0924;
      else
        return 0.1484;
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

  // Get the photon ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  event.getByToken(m_phoLooseIdMapToken, loose_id_decisions);
  event.getByToken(m_phoMediumIdMapToken, medium_id_decisions);
  event.getByToken(m_phoTightIdMapToken, tight_id_decisions);

  reset();
  m_size = 0;

  //std::cout<<"write info photon before loop"<<std::endl;
  for (unsigned int i = 0; i < p_photons.size(); ++i)
  {
    //std::cout<<"photon #"<<i<<std::endl;
    pat::PhotonRef photonRef(photonHandle, i);    

    // Look up and save the ID decisions
    m_pho_passLooseId[i] = (*loose_id_decisions)[photonRef];
    m_pho_passMediumId[i] = (*medium_id_decisions)[photonRef];
    m_pho_passTightId[i] = (*tight_id_decisions)[photonRef];

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

  edm::Handle<double> rhos;
  event.getByToken(m_rhoToken, rhos);
  double rho = *rhos;

  // All of the above can be replaced with
  const double chIso03 = part.chargedHadronIso();
  const double nhIso03 = part.neutralHadronIso();
  const double phIso03 = part.photonIso();


  // Now, isolations
  m_pho_chargedHadronsIsolation[index] = getCorrectedPFIsolation(chIso03, rho, photonRef->eta(), IsolationType::CHARGED_HADRONS);
  m_pho_neutralHadronsIsolation[index] = getCorrectedPFIsolation(nhIso03, rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS);
  m_pho_photonIsolation[index]         = getCorrectedPFIsolation(phIso03, rho, photonRef->eta(), IsolationType::PHOTONS);

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
    m_pho_chargedHadronsIsolation[i] = -1;
    m_pho_neutralHadronsIsolation[i] = -1;
    m_pho_photonIsolation[i] = -1;
    m_pho_passLooseId[i] = false;
    m_pho_passMediumId[i] = false;
    m_pho_passTightId[i] = false;
    m_pho_MCIndex[i] = -1;
  }
  
  if(m_pho_lorentzvector)
    m_pho_lorentzvector->Clear();
}


void PhotonExtractor::fillTree()
{
  m_tree_photon->Fill(); 
}


