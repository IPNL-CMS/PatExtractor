#include "../interface/MuonExtractor.h"

#include <DataFormats/VertexReco/interface/Vertex.h>

MuonExtractor::MuonExtractor(const std::string& name, const edm::InputTag& tag, const edm::InputTag& vertexTag, bool doTree)
  : BaseExtractor(name)
{
  m_OK = false;
  m_vertexTag = vertexTag;

  // Set everything to 0
 
  setPF((tag.label()).find("PFlow")); 

  m_muo_lorentzvector = new TClonesArray("TLorentzVector");
  m_scaleFactors.setWriteMode();

  reset();

  m_tag = tag;
  m_deltaR_cut = 0.5; // Maximum acceptable distance for MC matching

  // Tree definition

  if (doTree)
  {
    m_OK = true;

    m_tree_muon         = new TTree(m_name.c_str(), "PAT PF muon info"); 
    m_tree_muon->Branch("n_muons",  &m_size,  "n_muons/i");  
    m_tree_muon->Branch("muon_4vector","TClonesArray",&m_muo_lorentzvector, 1000, 0);
    m_tree_muon->Branch("muon_vx",  &m_muo_vx,   "muon_vx[n_muons]/F");  
    m_tree_muon->Branch("muon_vy",  &m_muo_vy,   "muon_vy[n_muons]/F");  
    m_tree_muon->Branch("muon_vz",  &m_muo_vz,   "muon_vz[n_muons]/F");  
    m_tree_muon->Branch("muon_charge", &m_muo_charge,  "muon_charge[n_muons]/I");
    m_tree_muon->Branch("muon_isGlobal", 	&m_muo_isGlobal,  "muon_isGlobal[n_muons]/I");
    m_tree_muon->Branch("muon_isTracker", &m_muo_isTracker, "muon_isTracker[n_muons]/I");
    m_tree_muon->Branch("muon_dB",        &m_muo_dB,        "muon_dB[n_muons]/F");
    m_tree_muon->Branch("muon_normChi2",  &m_muo_normChi2,  "muon_normChi2[n_muons]/F");
    m_tree_muon->Branch("muon_nValTrackerHits",&m_muo_nValTrackerHits,"muon_nValTrackerHits[n_muons]/I");
    m_tree_muon->Branch("muon_nValPixelHits",  &m_muo_nValPixelHits,"muon_nValPixelHits[n_muons]/I");
    m_tree_muon->Branch("muon_nMatches",       &m_muo_nMatches,"muon_nMatches[n_muons]/I");
    m_tree_muon->Branch("muon_trackIso",       &m_muo_trackIso,"muon_trackIso[n_muons]/F");
    m_tree_muon->Branch("muon_ecalIso",        &m_muo_ecalIso,"muon_ecalIso[n_muons]/F");
    m_tree_muon->Branch("muon_hcalIso",        &m_muo_hcalIso,"muon_hcalIso[n_muons]/F");
    m_tree_muon->Branch("muon_pfParticleIso",      &m_muo_pfParticleIso,"muon_pfParticleIso[n_muons]/F");
    m_tree_muon->Branch("muon_pfChargedHadronIso", &m_muo_pfChargedHadronIso,"muon_pfChargedHadronIso[n_muons]/F");
    m_tree_muon->Branch("muon_pfNeutralHadronIso", &m_muo_pfNeutralHadronIso,"muon_pfNeutralHadronIso[n_muons]/F");
    m_tree_muon->Branch("muon_pfPhotonIso",        &m_muo_pfPhotonIso,"muon_pfPhotonIso[n_muons]/F");
    m_tree_muon->Branch("muon_d0",      &m_muo_d0,"muon_d0[n_muons]/F");
    m_tree_muon->Branch("muon_d0error", &m_muo_d0error,"muon_d0error[n_muons]/F");
    m_tree_muon->Branch("muon_mcParticleIndex",&m_muo_MCIndex,"muon_mcParticleIndex[n_muons]/I");

    m_tree_muon->Branch("muon_nMatchedStations",              &m_muo_nMatchedStations,             "muon_nMatches[n_muons]/I");
    m_tree_muon->Branch("muon_trackerLayersWithMeasurement",  &m_muo_trackerLayersWithMeasurement, "muon_trackerLayersWithMeasurement[n_muons]/F");
    m_tree_muon->Branch("muon_dZ",                            &m_muo_dZ,                           "muon_dZ[n_muons]/F");
    m_tree_muon->Branch("muon_pixelLayerWithMeasurement",     &m_muo_pixelLayerWithMeasurement,    "muon_pixelLayerWithMeasurement[n_muons]/F");
    m_tree_muon->Branch("muon_globalTrackNumberOfValidHits",  &m_muo_globalTrackNumberOfValidHits, "muon_globalTrackNumberOfValidHits[n_muons]/F");

    m_tree_muon->Branch("muon_relIsolation",                   &m_muo_relIsolation, "muon_relIsolation[n_muons]/F");
    m_tree_muon->Branch("muon_deltaBetaCorrectedRelIsolation", &m_muo_deltaBetaCorrectedRelIsolation, "muon_deltaBetaCorrectedRelIsolation[n_muons]/F");

    m_tree_muon->Branch("muon_scaleFactor", &m_scaleFactors.getBackingArray());
  }

}

MuonExtractor::MuonExtractor(const std::string& name, TFile *a_file)
  : BaseExtractor(name)
{
  m_file = a_file;
  std::cout << "MuonExtractor objet is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_muon = dynamic_cast<TTree*>(a_file->Get(m_name.c_str()));

  if (!m_tree_muon)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_muo_lorentzvector = new TClonesArray("TLorentzVector");

  // Branches definition

  if (m_tree_muon->FindBranch("n_muons"))
    m_tree_muon->SetBranchAddress("n_muons",  &m_size);
  if (m_tree_muon->FindBranch("muon_4vector"))
    m_tree_muon->SetBranchAddress("muon_4vector",&m_muo_lorentzvector);
  if (m_tree_muon->FindBranch("muon_vx"))
    m_tree_muon->SetBranchAddress("muon_vx",  &m_muo_vx);
  if (m_tree_muon->FindBranch("muon_vy"))
    m_tree_muon->SetBranchAddress("muon_vy",  &m_muo_vy);
  if (m_tree_muon->FindBranch("muon_vz"))
    m_tree_muon->SetBranchAddress("muon_vz",  &m_muo_vz);
  if (m_tree_muon->FindBranch("muon_charge"))
  m_tree_muon->SetBranchAddress("muon_charge", &m_muo_charge);
  if (m_tree_muon->FindBranch("muon_isGlobal"))
    m_tree_muon->SetBranchAddress("muon_isGlobal", 	&m_muo_isGlobal);
  if (m_tree_muon->FindBranch("muon_isTracker"))
    m_tree_muon->SetBranchAddress("muon_isTracker", &m_muo_isTracker);
  if (m_tree_muon->FindBranch("muon_dB"))
    m_tree_muon->SetBranchAddress("muon_dB",        &m_muo_dB);
  if (m_tree_muon->FindBranch("muon_normChi2"))
    m_tree_muon->SetBranchAddress("muon_normChi2",  &m_muo_normChi2);
  if (m_tree_muon->FindBranch("muon_nValTrackerHits"))
    m_tree_muon->SetBranchAddress("muon_nValTrackerHits",&m_muo_nValTrackerHits);
  if (m_tree_muon->FindBranch("muon_nValPixelHits"))
    m_tree_muon->SetBranchAddress("muon_nValPixelHits",  &m_muo_nValPixelHits);
  if (m_tree_muon->FindBranch("muon_nMatches"))
    m_tree_muon->SetBranchAddress("muon_nMatches",       &m_muo_nMatches);
  if (m_tree_muon->FindBranch("muon_trackIso"))
    m_tree_muon->SetBranchAddress("muon_trackIso",       &m_muo_trackIso);
  if (m_tree_muon->FindBranch("muon_ecalIso"))
    m_tree_muon->SetBranchAddress("muon_ecalIso",        &m_muo_ecalIso);
  if (m_tree_muon->FindBranch("muon_hcalIso"))
    m_tree_muon->SetBranchAddress("muon_hcalIso",        &m_muo_hcalIso);
  if (m_tree_muon->FindBranch("muon_pfParticleIso"))
    m_tree_muon->SetBranchAddress("muon_pfParticleIso",      &m_muo_pfParticleIso);
  if (m_tree_muon->FindBranch("muon_pfChargedHadronIso"))
    m_tree_muon->SetBranchAddress("muon_pfChargedHadronIso", &m_muo_pfChargedHadronIso);
  if (m_tree_muon->FindBranch("muon_pfNeutralHadronIso"))
    m_tree_muon->SetBranchAddress("muon_pfNeutralHadronIso", &m_muo_pfNeutralHadronIso);
  if (m_tree_muon->FindBranch("muon_pfPhotonIso"))
    m_tree_muon->SetBranchAddress("muon_pfPhotonIso",        &m_muo_pfPhotonIso);
  if (m_tree_muon->FindBranch("muon_d0"))
    m_tree_muon->SetBranchAddress("muon_d0",      &m_muo_d0);
  if (m_tree_muon->FindBranch("muon_d0error"))
    m_tree_muon->SetBranchAddress("muon_d0error", &m_muo_d0error);
  if (m_tree_muon->FindBranch("muon_mcParticleIndex"))
    m_tree_muon->SetBranchAddress("muon_mcParticleIndex",&m_muo_MCIndex);  

  if (m_tree_muon->FindBranch("muon_nMatchedStations"))
    m_tree_muon->SetBranchAddress("muon_nMatchedStations", &m_muo_nMatchedStations);

  if (m_tree_muon->FindBranch("muon_trackerLayersWithMeasurement"))
    m_tree_muon->SetBranchAddress("muon_trackerLayersWithMeasurement", &m_muo_trackerLayersWithMeasurement);

  if (m_tree_muon->FindBranch("muon_dZ"))
    m_tree_muon->SetBranchAddress("muon_dZ", &m_muo_dZ);

  if (m_tree_muon->FindBranch("muon_pixelLayerWithMeasurement"))
    m_tree_muon->SetBranchAddress("muon_pixelLayerWithMeasurement", &m_muo_pixelLayerWithMeasurement);

  if (m_tree_muon->FindBranch("muon_globalTrackNumberOfValidHits"))
    m_tree_muon->SetBranchAddress("muon_globalTrackNumberOfValidHits", &m_muo_globalTrackNumberOfValidHits);

  if (m_tree_muon->FindBranch("muon_relIsolation"))
    m_tree_muon->SetBranchAddress("muon_relIsolation", &m_muo_relIsolation);

  if (m_tree_muon->FindBranch("muon_deltaBetaCorrectedRelIsolation"))
    m_tree_muon->SetBranchAddress("muon_deltaBetaCorrectedRelIsolation", &m_muo_deltaBetaCorrectedRelIsolation);

  if (m_tree_muon->FindBranch("muon_scaleFactor"))
    m_tree_muon->SetBranchAddress("muon_scaleFactor", &m_scaleFactors.getBackingArray());
}

MuonExtractor::~MuonExtractor()
{
  delete m_muo_lorentzvector;
}

//
// Method getting the info from an input file
//

void MuonExtractor::getInfo(int ievt) 
{
  m_tree_muon->GetEntry(ievt); 
}

void MuonExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Muon& part, int index) 
{
  if (index>=m_muons_MAX) return;

  edm::Handle<std::vector<reco::Vertex>> pvHandle;
  event.getByLabel(m_vertexTag, pvHandle);

  new((*m_muo_lorentzvector)[index]) TLorentzVector(part.px(),part.py(),part.pz(),part.energy());
  m_muo_vx[index]                             = part.vx();
  m_muo_vy[index]                             = part.vy();
  m_muo_vz[index]                             = part.vz();
  m_muo_isGlobal[index]                       = part.isGlobalMuon();
  m_muo_isTracker[index]                      = part.isTrackerMuon();
  m_muo_charge[index]                         = part.charge();
  m_muo_nMatchedStations[index]               = part.numberOfMatchedStations();
  m_muo_dZ[index]                             = fabs(part.muonBestTrack()->dz(pvHandle->at(0).position()));
  m_muo_trackerLayersWithMeasurement[index]   = part.track()->hitPattern().trackerLayersWithMeasurement();
  
  if (part.outerTrack().isNonnull())
  {
    m_muo_dB[index]              = part.dB();
    m_muo_nMatches[index]        = part.numberOfMatches();
  }
  
  if (part.innerTrack().isNonnull())
  {
    m_muo_dB[index]                        = part.dB();
    m_muo_pixelLayerWithMeasurement[index] = part.innerTrack()->hitPattern().pixelLayersWithMeasurement();
    m_muo_nValPixelHits[index]             = part.innerTrack()->hitPattern().numberOfValidPixelHits();
    m_muo_nValTrackerHits[index]           = part.innerTrack()->hitPattern().numberOfValidTrackerHits();
    m_muo_d0[index]                        = part.innerTrack()->d0(); 
    m_muo_d0error[index]                   = part.innerTrack()->d0Error();
  }
  
  if (part.globalTrack().isNonnull())
  {
    m_muo_nValTrackerHits[index]              = part.numberOfValidHits();
    m_muo_normChi2[index]                     = part.normChi2();
    m_muo_globalTrackNumberOfValidHits[index] = part.globalTrack()->hitPattern().numberOfValidMuonHits();
  }

  m_muo_trackIso[index]        = part.trackIso();
  m_muo_ecalIso[index]         = part.ecalIso();
  m_muo_hcalIso[index]         = part.hcalIso();
  
  if (m_isPF)
  {
    m_muo_pfParticleIso[index]      = part.particleIso();
    m_muo_pfChargedHadronIso[index] = part.chargedHadronIso();
    m_muo_pfNeutralHadronIso[index] = part.neutralHadronIso();
    m_muo_pfPhotonIso[index]        = part.photonIso();

    m_muo_relIsolation[index] = (part.chargedHadronIso() + part.neutralHadronIso() + part.photonIso()) / part.pt();
    m_muo_deltaBetaCorrectedRelIsolation[index] = (part.chargedHadronIso() + std::max((part.neutralHadronIso() + part.photonIso()) - 0.5 * part.puChargedHadronIso(), 0.0)) / part.pt();
  } 

  if (m_isMC)
    m_scaleFactors.push_back(m_scaleFactorService->getMuonScaleFactor(part.pt(), part.eta()));
}


// Method initializing everything (to do for each event)

void MuonExtractor::reset()
{
  m_size = 0;
  m_scaleFactors.clear();

  for (int i=0;i<m_muons_MAX;++i) 
  {
    m_muo_vx[i] = 0.;
    m_muo_vy[i] = 0.;
    m_muo_vz[i] = 0.;
    m_muo_isGlobal[i] = 0;
    m_muo_isTracker[i] = 0;
    m_muo_charge[i] = 0;
    m_muo_dB[i] = 0.;
    m_muo_normChi2[i] = 0.;
    m_muo_nValTrackerHits[i] = 0;
    m_muo_nValPixelHits[i] = 0;
    m_muo_nMatches[i] = 0;
    m_muo_trackIso[i] = 0.;
    m_muo_ecalIso[i] = 0.;
    m_muo_hcalIso[i] = 0.;
    m_muo_pfParticleIso[i] =0.;
    m_muo_pfChargedHadronIso[i]=0.;
    m_muo_pfNeutralHadronIso[i]=0.;
    m_muo_pfPhotonIso[i]=0.;
    m_muo_d0[i]=0.;
    m_muo_d0error[i]=0.;
    m_muo_MCIndex[i]=-1;

    m_muo_nMatchedStations[i] = -1;
    m_muo_trackerLayersWithMeasurement[i] = -1;
    m_muo_dZ[i] = -1;
    m_muo_pixelLayerWithMeasurement[i] = -1;
    m_muo_globalTrackNumberOfValidHits[i] = -1;

    m_muo_relIsolation[i] = -1;
    m_muo_deltaBetaCorrectedRelIsolation[i] = -1;
  }
  m_muo_lorentzvector->Clear();
}


void MuonExtractor::fillTree()
{
  m_tree_muon->Fill(); 
}
