#include "../interface/PFpartExtractor.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/MuonReco/interface/Muon.h"

using namespace std;

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
    m_tree_pfpart->Branch("pf_4vector","TClonesArray",&m_pf_lorentzvector, 1000, 0);
    m_tree_pfpart->Branch("pf_vx",     &m_pf_vx,     "pf_vx[n_pf]/F");  
    m_tree_pfpart->Branch("pf_vy",     &m_pf_vy,     "pf_vy[n_pf]/F");  
    m_tree_pfpart->Branch("pf_vz",     &m_pf_vz,     "pf_vz[n_pf]/F");  
    m_tree_pfpart->Branch("pf_charge", &m_pf_charge, "pf_charge[n_pf]/I");
    m_tree_pfpart->Branch("pf_pdgid",  &m_pf_pdgid,  "pf_charge[n_pf]/I");
    m_tree_pfpart->Branch("pf_trkLayer", &m_pf_trkLayer,  "pf_trkLayer[n_pf]/I");
    m_tree_pfpart->Branch("pf_pixLayer", &m_pf_pixLayer,  "pf_pixLayer[n_pf]/I");
    m_tree_pfpart->Branch("pf_trknormChi2", &m_pf_trknormChi2,  "pf_trknormChi2[n_pf]/F");
    m_tree_pfpart->Branch("pf_numberOfChambers", &m_pf_numberOfChambers,  "pf_numberOfChambers[n_pf]/I");
    m_tree_pfpart->Branch("pf_numberOfMatchedStations", &m_pf_numberOfMatchedStations,  "pf_numberOfMatchedStations[n_pf]/I");
    m_tree_pfpart->Branch("pf_isGlobalMuon", &m_pf_isGlobalMuon,  "pf_isGlobalMuon[n_pf]/O");
    m_tree_pfpart->Branch("pf_isTrackerMuon", &m_pf_isTrackerMuon,  "pf_isTrackerMuon[n_pf]/O");
    m_tree_pfpart->Branch("pf_isStandAloneMuon", &m_pf_isStandAloneMuon,  "pf_isStandAloneMuon[n_pf]/O");
    m_tree_pfpart->Branch("pf_isCaloMuon", &m_pf_isCaloMuon,  "pf_isCaloMuon[n_pf]/O");
    m_tree_pfpart->Branch("pf_isPFMuon", &m_pf_isPFMuon,  "pf_isPFMuon[n_pf]/O");
    m_tree_pfpart->Branch("pf_isRPCMuon", &m_pf_isRPCMuon,  "pf_isRPCMuon[n_pf]/O");

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

  m_tree_pfpart->SetBranchAddress("n_pf",     &m_pf_size);
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

// Dummy pure virtual function
void PFpartExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const reco::PFCandidate& part, int index) 
{

}

void PFpartExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor)
{
  edm::Handle<reco::PFCandidateCollection>  pfHandle;
  event.getByLabel(m_tag, pfHandle);
  reco::PFCandidateCollection pfs = *pfHandle;
  if ( !pfHandle.isValid() ) {
    std::cout << "PFpartExtractor::writeInfo(): pfHandle is not valid..." << std::endl;
    return;
  }

  // Select good PF muons :
  //-----------------------

  reset();

  m_size = 0;

  std::vector<const reco::PFCandidate*> myPFparts;
  std::vector<int> indpfs;

  for ( unsigned int i = 0; i < pfs.size(); ++i ) {

    if ( abs(pfs[i].pdgId()) != 13  ) continue;
    if ( pfs[i].pt()         <   4. ) continue;

    myPFparts.push_back(&pfs[i]);

    if (m_pf_size>=m_pfpart_MAX) continue;

    new((*m_pf_lorentzvector)[m_pf_size]) TLorentzVector(pfs[i].px(),pfs[i].py(),pfs[i].pz(),pfs[i].energy());
    if ( !pfs[i].trackRef() ) {
      m_pf_vx[m_pf_size]              = 0;
      m_pf_vy[m_pf_size]              = 0;
      m_pf_vz[m_pf_size]              = 0;
      m_pf_trkLayer[m_pf_size]        = 0;
      m_pf_pixLayer[m_pf_size]        = 0;
    } else {
      m_pf_vx[m_pf_size]              = pfs[i].trackRef()->vx();
      m_pf_vy[m_pf_size]              = pfs[i].trackRef()->vy();
      m_pf_vz[m_pf_size]              = pfs[i].trackRef()->vz();
      m_pf_trkLayer[m_pf_size]        = pfs[i].trackRef()->hitPattern().trackerLayersWithMeasurement();
      m_pf_pixLayer[m_pf_size]        = pfs[i].trackRef()->hitPattern().pixelLayersWithMeasurement();
      m_pf_trknormChi2[m_pf_size]     = pfs[i].trackRef()->normalizedChi2();
    }

    if ( !pfs[i].muonRef() ) {
      cout << "PF muonRef is invalid ... " << endl;
    } else {
      //      m_pf_numberOfChambers[m_pf_size]        = pfs[i].muonRef()->numberOfChambers();
      //      m_pf_numberOfMatchedStations[m_pf_size] = pfs[i].muonRef()->numberOfMatchedStations();
      //      m_pf_isGlobalMuon[m_pf_size]     =  pfs[i].muonRef()->isGlobalMuon();
      //      m_pf_isTrackerMuon[m_pf_size]    =  pfs[i].muonRef()->isTrackerMuon();
      //      m_pf_isStandAloneMuon[m_pf_size] =  pfs[i].muonRef()->isStandAloneMuon();
      //      m_pf_isCaloMuon[m_pf_size]       =  pfs[i].muonRef()->isCaloMuon();
      //      m_pf_isPFMuon[m_pf_size]         =  pfs[i].muonRef()->isPFMuon();
      //      m_pf_isRPCMuon[m_pf_size]        =  pfs[i].muonRef()->isRPCMuon();
    }

    m_pf_charge[m_pf_size]          = pfs[i].charge();
    m_pf_pdgid[m_pf_size]           = pfs[i].pdgId();


    indpfs.push_back(m_pf_size);
    ++m_pf_size;
  }

  fillTree();
}


// Method initializing everything (to do for each event)

void PFpartExtractor::reset()
{

  // PF particle tree

  m_size = 0;
  m_pf_size = 0;

  for (int i=0;i<m_pfpart_MAX;++i) {
    m_pf_vx[i]         = 0.;
    m_pf_vy[i]         = 0.;
    m_pf_vz[i]         = 0.;
    m_pf_charge[i]     = 0;
    m_pf_pdgid[i]      = 0.;
    m_pf_trkLayer[i]   = 0;
    m_pf_pixLayer[i]   = 0;
    m_pf_trknormChi2[i]= 0.;
    m_pf_numberOfChambers[i]        = 0;
    m_pf_numberOfMatchedStations[i] = 0;
    m_pf_isGlobalMuon[i] = false;
    m_pf_isTrackerMuon[i] = false;
    m_pf_isStandAloneMuon[i] = false;
    m_pf_isCaloMuon[i] = false;
    m_pf_isPFMuon[i] = false;
    m_pf_isRPCMuon[i] = false;
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

