#include "../interface/MuonExtractor.h"


MuonExtractor::MuonExtractor(bool doTree, edm::InputTag tag)
{
  m_isPF_muon=true; // By default we use PF muons
  m_OK = false;

  // Set everything to 0
 
  MuonExtractor::setPF((tag.label()).find("PFlow")); 
  m_muo_lorentzvector = new TClonesArray("TLorentzVector");
  MuonExtractor::reset();

  m_tag = tag;
  m_deltaR_cut = 0.5; // Maximum acceptable distance for MC matching


  // Tree definition

  if (doTree)
  {
    m_OK = true;

    m_tree_muon         = new TTree("muon_PF","PAT PF muon info"); 
    m_tree_muon->Branch("n_muons",  &m_n_muons,  "n_muons/I");  
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
  }

}

MuonExtractor::MuonExtractor(TFile *a_file)
{
  std::cout << "MuonExtractor objet is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_muon = dynamic_cast<TTree*>(a_file->Get("muon_PF"));

  if (!m_tree_muon)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_muo_lorentzvector = new TClonesArray("TLorentzVector");


  // Branches definition

  if (m_tree_muon->FindBranch("n_muons"))
    m_tree_muon->SetBranchAddress("n_muons",  &m_n_muons);
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
}

MuonExtractor::~MuonExtractor()
{}



//
// Method filling the main particle tree
//

void MuonExtractor::writeInfo(const edm::Event *event,MCExtractor* m_MC, bool doMC) 
{
  edm::Handle< edm::View<pat::Muon> >  muonHandle;
  event->getByLabel(m_tag, muonHandle);
  edm::View<pat::Muon> p_muons = *muonHandle;

  MuonExtractor::reset();
  MuonExtractor::fillSize(static_cast<int>(p_muons.size()));
  
  if (MuonExtractor::getSize())
  {
    for(int i=0; i<MuonExtractor::getSize(); ++i) 
    {
      MuonExtractor::writeInfo(&p_muons.at(i),i);
      
      if (doMC)
      {
	int idx_min = MuonExtractor::getMatch(&p_muons.at(i),m_MC); 
	m_muo_MCIndex[i] = idx_min;
      }
    }
  }

  MuonExtractor::fillTree();
}

//
// Method getting the info from an input file
//

void MuonExtractor::getInfo(int ievt) 
{
  m_tree_muon->GetEntry(ievt); 
}

void MuonExtractor::writeInfo(const pat::Muon *part, int index) 
{
  if (index>=m_muons_MAX) return;

  new((*m_muo_lorentzvector)[index]) TLorentzVector(part->px(),part->py(),part->pz(),part->energy());
  m_muo_vx[index]              = part->vx();
  m_muo_vy[index]              = part->vy();
  m_muo_vz[index]              = part->vz();
  m_muo_isGlobal[index]        = part->isGlobalMuon();
  m_muo_isTracker[index]       = part->isTrackerMuon();
  m_muo_charge[index]          = part->charge();
  
  if (part->outerTrack().isNonnull())
  {
    m_muo_dB[index]              = part->dB();
    m_muo_nMatches[index]        = part->numberOfMatches();
  }
  
  if (part->innerTrack().isNonnull())
  {
    m_muo_dB[index]                = part->dB();
    m_muo_nValPixelHits[index]     = part->innerTrack()->hitPattern().pixelLayersWithMeasurement();
    m_muo_nValTrackerHits[index]   = part->innerTrack()->hitPattern().numberOfValidTrackerHits();
    m_muo_d0[index]                = part->innerTrack()->d0(); 
    m_muo_d0error[index]           = part->innerTrack()->d0Error();
  }
  
  if (part->globalTrack().isNonnull())
  {
    m_muo_nValTrackerHits[index]   = part->numberOfValidHits();
    m_muo_normChi2[index]          = part->normChi2();
  }

  m_muo_trackIso[index]        = part->trackIso();
  m_muo_ecalIso[index]         = part->ecalIso();
  m_muo_hcalIso[index]         = part->hcalIso();
  
  if (m_isPF_muon)
  {
    m_muo_pfParticleIso[index]      = part->particleIso();
    m_muo_pfChargedHadronIso[index] = part->chargedHadronIso();
    m_muo_pfNeutralHadronIso[index] = part->neutralHadronIso();
    m_muo_pfPhotonIso[index]        = part->photonIso();
  } 
}


// Method initializing everything (to do for each event)

void MuonExtractor::reset()
{
  m_n_muons = 0;

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
  }
  m_muo_lorentzvector->Clear();
}


void MuonExtractor::fillTree()
{
  m_tree_muon->Fill(); 
}
 
void MuonExtractor::fillSize(int size)
{
  m_n_muons=size;
}

int  MuonExtractor::getSize()
{
  return m_n_muons;
}

void MuonExtractor::setPF(bool isPF)
{
  m_isPF_muon=isPF;
}

int MuonExtractor::getMatch(const pat::Muon *part, MCExtractor* m_MC)
{
  float deltaR_min = 1e6;
  int idx_min    = -1;
  float deltaR;    
  float deltaP;
  TLorentzVector TL_genPart;
  TLorentzVector TL_muon;

  for(int mcPart_i=0; mcPart_i<m_MC->getSize(); ++mcPart_i) 
  {
    if (m_MC->getStatus(mcPart_i)!=3) continue;
    if (fabs(m_MC->getType(mcPart_i))!=13) continue;
    
    TL_genPart.SetPxPyPzE(m_MC->getPx(mcPart_i),m_MC->getPy(mcPart_i),m_MC->getPz(mcPart_i),m_MC->getE(mcPart_i));
    TL_muon.SetPxPyPzE(part->px(),part->py(),part->pz(),part->energy());
    
    if(TL_genPart.Pt())
    {
      deltaR = TL_genPart.DeltaR(TL_muon);
      deltaP = fabs(TL_genPart.Pt()-TL_muon.Pt())/(TL_muon.Pt());
      if (deltaP>0.5) continue; //min DPrel for muons 0.5
      if (deltaR>0.5) continue; //min DR for muons 0.5
      //      if((m_MC->getType(mcPart_i)*part->charge())<0.) continue;
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
