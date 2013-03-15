#include "../interface/TrackExtractor.h"


TrackExtractor::TrackExtractor(const std::string& name, const edm::InputTag& tag, bool doTree)
  : BaseExtractor(name)
{
  m_tag = tag;

  // Set everything to 0
  m_OK = false;
  reset();

  // Tree definition

  if (doTree)
  {
    m_OK = true;
    m_tree_track         = new TTree(name.c_str(), "General tracks info");     
    m_tree_track->Branch("n_tracks",  &m_size,  "n_tracks/i");  
    m_tree_track->Branch("track_px",  &m_trk_px,   "track_px[n_tracks]/F");  
    m_tree_track->Branch("track_py",  &m_trk_py,   "track_py[n_tracks]/F");  
    m_tree_track->Branch("track_pz",  &m_trk_pz,   "track_pz[n_tracks]/F"); 
    m_tree_track->Branch("track_vx",  &m_trk_vx,   "track_vx[n_tracks]/F");  
    m_tree_track->Branch("track_vy",  &m_trk_vy,   "track_vy[n_tracks]/F");  
    m_tree_track->Branch("track_vz",  &m_trk_vz,   "track_vz[n_tracks]/F");  
    m_tree_track->Branch("track_eta", &m_trk_eta,  "track_eta[n_tracks]/F");  
    m_tree_track->Branch("track_phi", &m_trk_phi,  "track_phi[n_tracks]/F");
    m_tree_track->Branch("track_charge", &m_trk_charge,  "track_charge[n_tracks]/I");
    m_tree_track->Branch("track_d0",        &m_trk_d0,        "track_d0[n_tracks]/F");
    m_tree_track->Branch("track_normChi2",  &m_trk_normChi2,  "track_normChi2[n_tracks]/F");
    m_tree_track->Branch("track_nValidHits",&m_trk_nValidHits,"track_nValidHits[n_tracks]/I");
  }
}

TrackExtractor::TrackExtractor(const std::string& name, TFile *a_file)
  :BaseExtractor(name)
{
  std::cout << "TrackExtractor objet is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_track = dynamic_cast<TTree*>(a_file->Get(m_name.c_str()));

  if (!m_tree_track)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_tree_track->SetBranchAddress("n_tracks",  &m_size);
  m_tree_track->SetBranchAddress("track_px",  &m_trk_px);
  m_tree_track->SetBranchAddress("track_py",  &m_trk_py);
  m_tree_track->SetBranchAddress("track_pz",  &m_trk_pz);
  m_tree_track->SetBranchAddress("track_vx",  &m_trk_vx);
  m_tree_track->SetBranchAddress("track_vy",  &m_trk_vy);
  m_tree_track->SetBranchAddress("track_vz",  &m_trk_vz);
  m_tree_track->SetBranchAddress("track_eta", &m_trk_eta);
  m_tree_track->SetBranchAddress("track_phi", &m_trk_phi);
  m_tree_track->SetBranchAddress("track_charge", &m_trk_charge);
  m_tree_track->SetBranchAddress("track_d0",        &m_trk_d0);
  m_tree_track->SetBranchAddress("track_normChi2",  &m_trk_normChi2);
  m_tree_track->SetBranchAddress("track_nValidHits",&m_trk_nValidHits);

}

TrackExtractor::~TrackExtractor()
{}



//
// Method filling the main particle tree
//

void TrackExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const reco::Track& part, int index) 
{
  if (index>=m_tracks_MAX) return;

  m_trk_px[index]              = part.px();
  m_trk_py[index]              = part.py();
  m_trk_pz[index]              = part.pz();
  m_trk_vx[index]              = part.vx();
  m_trk_vy[index]              = part.vy();
  m_trk_vz[index]              = part.vz();
  m_trk_eta[index]             = part.eta();
  m_trk_phi[index]             = part.phi();
  m_trk_charge[index]          = part.charge();
  m_trk_d0[index]              = part.d0();
  m_trk_nValidHits[index]      = part.numberOfValidHits();
  m_trk_normChi2[index]        = part.normalizedChi2();
}


// Method initializing everything (to do for each event)

void TrackExtractor::reset()
{
  m_size = 0;

  for (int i=0;i<m_tracks_MAX;++i) 
  {
    m_trk_px[i]         = 0.;
    m_trk_py[i]         = 0.;
    m_trk_pz[i]         = 0.;
    m_trk_vx[i]         = 0.;
    m_trk_vy[i]         = 0.;
    m_trk_vz[i]         = 0.;
    m_trk_eta[i]        = 0.;
    m_trk_phi[i]        = 0.;
    m_trk_charge[i]     = 0;
    m_trk_d0[i]         = 0.;
    m_trk_normChi2[i]   = 0.;
    m_trk_nValidHits[i] = 0;
  }
}


void TrackExtractor::fillTree()
{
  m_tree_track->Fill(); 
}
 
//
// Method getting the info from an input file
//

void TrackExtractor::getInfo(int ievt) 
{
  m_tree_track->GetEntry(ievt); 
}
