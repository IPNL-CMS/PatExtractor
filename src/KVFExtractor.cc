#include "../interface/KVFExtractor.h"
#include "../interface/JECReader.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>

#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
// #include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#define DEBUG false

using namespace std;

//----------------------------------------------------------------------
// functions from : RecoVertex/KinematicFit/plugins/KineExample.cc

void printOut(const RefCountedKinematicVertex& myVertex)
{
  if (myVertex->vertexIsValid()) {
    cout << "Decay vertex: " << myVertex->position() <<myVertex->chiSquared()<< " "<<myVertex->degreesOfFreedom()<<endl;
  } else cout << "Decay vertex Not valid\n";
}

void printOut(const RefCountedKinematicParticle& myParticle)
{
  cout << "Particle: \n";
  //accessing the reconstructed Bs meson parameters:
  //SK: uncomment if needed  AlgebraicVector7 bs_par = myParticle->currentState().kinematicParameters().vector();

  //and their joint covariance matrix:
  //SK:uncomment if needed  AlgebraicSymMatrix77 bs_er = myParticle->currentState().kinematicParametersError().matrix();
  cout << "Momentum at vertex: " << myParticle->currentState().globalMomentum ()<<endl;
  cout << "Parameters at vertex: " << myParticle->currentState().kinematicParameters().vector()<<endl;
}

void printOut(const RefCountedKinematicTree& myTree)
{
  if (!myTree->isValid()) {
    cout <<"Tree is invalid. Fit failed.\n";
    return;
  }

  //accessing the tree components, move pointer to top
  myTree->movePointerToTheTop();

  //We are now at the top of the decay tree getting the jpsi reconstructed KinematicPartlcle
  RefCountedKinematicParticle jpsi = myTree->currentParticle();
  printOut(jpsi);

  // The jpsi decay vertex
  RefCountedKinematicVertex b_dec_vertex = myTree->currentDecayVertex();
  printOut(b_dec_vertex);

  // Get all the children of Jpsi:
  //In this way, the pointer is not moved
  vector< RefCountedKinematicParticle > jpsi_children = myTree->finalStateParticles();

  for (unsigned int i=0;i< jpsi_children.size();++i) {
    printOut(jpsi_children[i]);
  }

  //Now navigating down the tree , pointer is moved:
  bool child = myTree->movePointerToTheFirstChild();

  if(child) while (myTree->movePointerToTheNextChild()) {
    RefCountedKinematicParticle aChild = myTree->currentParticle();
    printOut(aChild);
  }
}

//-------------------------------------------------------------------------

  KVFExtractor::KVFExtractor(const std::string& name_jpsi, const std::string& name_d0, const std::string& name_jet, std::shared_ptr<ScaleFactorService> sf, const edm::ParameterSet& config)
: BaseExtractor(name_jet, sf)
{
  if (! config.exists(name_jpsi))
    throw edm::Exception(edm::errors::ConfigFileReadError) << "No edm::ParameterSet named " << name_jpsi << " found";
  if (! config.exists(name_d0))
    throw edm::Exception(edm::errors::ConfigFileReadError) << "No edm::ParameterSet named " << name_d0 << " found";
  if (! config.exists(name_jet))
    throw edm::Exception(edm::errors::ConfigFileReadError) << "No edm::ParameterSet named " << name_jet << " found";

  const edm::ParameterSet& jpsiConfig = config.getParameter<edm::ParameterSet>(name_jpsi);

  m_muJpsiMinPt = jpsiConfig.getUntrackedParameter<double>("muJpsiMinPt");
  m_jpsiMassMin = jpsiConfig.getUntrackedParameter<double>("jpsiMassMin");
  m_jpsiMassMax = jpsiConfig.getUntrackedParameter<double>("jpsiMassMax");
  m_vertexTag   = jpsiConfig.getUntrackedParameter<edm::InputTag>("vtxTag");

  const edm::ParameterSet& d0Config = config.getParameter<edm::ParameterSet>(name_d0);

  m_nTrD0Max = d0Config.getUntrackedParameter<unsigned int>("nTrD0Max");
  m_trSumMinPt = d0Config.getUntrackedParameter<double>("trSumMinPt");
  m_trUnfoldMinPt = d0Config.getUntrackedParameter<double>("trUnfoldMinPt");

  const edm::ParameterSet& jetConfig = config.getParameter<edm::ParameterSet>(name_jet);

  m_tag = jetConfig.getParameter<edm::InputTag>("input");

  mCorrectJets = jetConfig.getUntrackedParameter<bool>("redoJetCorrection", false);
  mUseGlobalTagForJEC = jetConfig.getUntrackedParameter<bool>("useGlobalTagForJEC", true);

  if (!mUseGlobalTagForJEC) {
    mJecPayload =  jetConfig.getUntrackedParameter<std::string>("jecPayload");
    mJecJetAlgo =  jetConfig.getUntrackedParameter<std::string>("jecJetAlgo");
    if (mJecPayload.length() > 0) {
      mJecPayload = edm::FileInPath(mJecPayload).fullPath();
    } else {
      std::cout << "WARNING! No JecPayload file found. Use the global tag instead for JEC" << std::endl;
      mUseGlobalTagForJEC = true;
    }
  }

  mTxtCorrector = nullptr;

  if (mCorrectJets)
    mJetCorrectorLabel = jetConfig.getParameter<std::string>("jetCorrectorLabel");
  mDoJER = jetConfig.getUntrackedParameter<bool>("doJER", true);
  mJERSign = 0;
  if (mDoJER) {
    mJERSign = jetConfig.getUntrackedParameter<int>("jerSign", 0);
  }
  mDoLooseJetID = jetConfig.getUntrackedParameter<bool>("doLooseJetID", true);

  mJESSign = jetConfig.getUntrackedParameter<int>("jesSign", 0);

  if (mJERSign != 0 && mJERSign != -1 && mJERSign != 1)
    throw edm::Exception(edm::errors::LogicError) << "jerSign must be 0 for nominal correction, -1 for 1-sigma down correction, or 1 for 1-sigma up correction";

  if (mJESSign != 0 && mJESSign != -1 && mJESSign != 1)
    throw edm::Exception(edm::errors::LogicError) << "jesSign must be 0 for nominal correction, -1 for 1-sigma down correction, or 1 for 1-sigma up correction";

  mJecFilename = jetConfig.getUntrackedParameter<std::string>("jes_uncertainties_file", "");
  if (mJESSign != 0 && mJecFilename.length() > 0) {
    mJecFilename = edm::FileInPath(mJecFilename).fullPath();
  }

#if DEBUG
  std::cout << "##########" << std::endl;
  std::cout << "KVF extractor summary" << std::endl;
  std::cout << "Redo JEC: " << mCorrectJets << ((mCorrectJets) ? ("; " + mJetCorrectorLabel) : "") << std::endl;
  std::cout << "Do JER: " << mDoJER << "; JER sign: " << mJERSign << std::endl;
#endif

  // Set everything to 0
  const auto& sfWorkingPoints = m_scaleFactorService->getMuonScaleFactorWorkingPoints();

  m_jpsi_jet_lorentzvector    = new TClonesArray("TLorentzVector");
  m_jpsi_jet_scaleFactors.setWriteMode();
  m_jpsipf_lorentzvector      = new TClonesArray("TLorentzVector");
  m_jpsikvf_lorentzvector     = new TClonesArray("TLorentzVector");
  m_jpsikvf_mu1_lorentzvector = new TClonesArray("TLorentzVector");
  for (auto& it: sfWorkingPoints) {
    std::string name = "muon_scaleFactor_" + ScaleFactorService::workingPointToString(it.first) + "eff_" + ScaleFactorService::workingPointToString(it.second) + "iso";
    m_jpsikvf_mu1_muon_scaleFactors[name] = ScaleFactorCollection();
    m_jpsikvf_mu1_muon_scaleFactors[name].setWriteMode();
  }
  m_jpsikvf_mu2_lorentzvector = new TClonesArray("TLorentzVector");
  for (auto& it: sfWorkingPoints) {
    std::string name = "muon_scaleFactor_" + ScaleFactorService::workingPointToString(it.first) + "eff_" + ScaleFactorService::workingPointToString(it.second) + "iso";
    m_jpsikvf_mu2_muon_scaleFactors[name] = ScaleFactorCollection();
    m_jpsikvf_mu2_muon_scaleFactors[name].setWriteMode();
  }

  m_mujet_jet_lorentzvector           = new TClonesArray("TLorentzVector");
  m_mujet_jet_scaleFactors.setWriteMode();
  m_mujet_nonisomuplus_lorentzvector  = new TClonesArray("TLorentzVector");
  for (auto& it: sfWorkingPoints) {
    std::string name = "muon_scaleFactor_" + ScaleFactorService::workingPointToString(it.first) + "eff_" + ScaleFactorService::workingPointToString(it.second) + "iso";
    m_mujet_nonisomuplus_scaleFactors[name] = ScaleFactorCollection();
    m_mujet_nonisomuplus_scaleFactors[name].setWriteMode();
  }
  m_mujet_nonisomuminus_lorentzvector = new TClonesArray("TLorentzVector");
  for (auto& it: sfWorkingPoints) {
    std::string name = "muon_scaleFactor_" + ScaleFactorService::workingPointToString(it.first) + "eff_" + ScaleFactorService::workingPointToString(it.second) + "iso";
    m_mujet_nonisomuminus_scaleFactors[name] = ScaleFactorCollection();
    m_mujet_nonisomuminus_scaleFactors[name].setWriteMode();
  }
  m_mujet_tr_lorentzvector            = new TClonesArray("TLorentzVector");
  m_mujet_d0pf_lorentzvector          = new TClonesArray("TLorentzVector");
  m_mujet_d0kvf_lorentzvector         = new TClonesArray("TLorentzVector");
  m_mujet_d0kvf_pion_lorentzvector    = new TClonesArray("TLorentzVector");
  m_mujet_d0kvf_kaon_lorentzvector    = new TClonesArray("TLorentzVector");

  reset();

  // Tree definition

  m_tree_jpsi = NULL;
  m_tree_jpsi = new TTree(name_jpsi.c_str(), "Jpsi info");  
  m_tree_jpsi->Branch("n_jpsi", &m_jpsi_size, "n_jpsi/I");
  m_tree_jpsi->Branch("jpsi_indjet", &m_jpsi_indjet, "jpsi_indjet[n_jpsi]/I");
  m_tree_jpsi->Branch("jpsi_jet_btag_CSV", &m_jpsi_jet_btag_CSV, "jpsi_jet_btag_CSV[n_jpsi]/F");
  m_tree_jpsi->Branch("jpsi_jet_4vector", "TClonesArray", &m_jpsi_jet_lorentzvector, 1000, 0);
  m_tree_jpsi->Branch("jpsi_jet_scaleFactor", &m_jpsi_jet_scaleFactors.getBackingArray());
  m_tree_jpsi->Branch("jpsi_indpf1", &m_jpsi_indpf1, "jpsi_indpf1[n_jpsi]/I");
  m_tree_jpsi->Branch("jpsi_indpf2", &m_jpsi_indpf2, "jpsi_indpf2[n_jpsi]/I");
  m_tree_jpsi->Branch("jpsipf_4vector", "TClonesArray", &m_jpsipf_lorentzvector, 1000, 0);
  m_tree_jpsi->Branch("jpsi_4vector", "TClonesArray", &m_jpsikvf_lorentzvector, 1000, 0);
  m_tree_jpsi->Branch("jpsi_mu1_4vector", "TClonesArray", &m_jpsikvf_mu1_lorentzvector, 1000, 0);
  for (auto& it: m_jpsikvf_mu1_muon_scaleFactors) {
    m_tree_jpsi->Branch(("jpsi_mu1_"+it.first).c_str(), & it.second.getBackingArray());
  }
  m_tree_jpsi->Branch("jpsi_mu2_4vector", "TClonesArray", &m_jpsikvf_mu2_lorentzvector, 1000, 0);
  for (auto& it: m_jpsikvf_mu2_muon_scaleFactors) {
    m_tree_jpsi->Branch(("jpsi_mu2_"+it.first).c_str(), & it.second.getBackingArray());
  }
  m_tree_jpsi->Branch("jpsi_vx", &m_jpsikvf_vx,	"jpsikvf_vx[n_jpsi]/F");  
  m_tree_jpsi->Branch("jpsi_vy", &m_jpsikvf_vy,	"jpsikvf_vy[n_jpsi]/F");  
  m_tree_jpsi->Branch("jpsi_vz", &m_jpsikvf_vz,	"jpsikvf_vz[n_jpsi]/F");
  m_tree_jpsi->Branch("jpsi_vtxvalid", &m_jpsikvf_vtxvalid, "jpsikvf_vtxvalid[n_jpsi]/O");
  m_tree_jpsi->Branch("jpsi_vtxchi2", &m_jpsikvf_vtxchi2,  "jpsikvf_vtxchi2[n_jpsi]/F");
  m_tree_jpsi->Branch("jpsi_ndf",	&m_jpsikvf_ndf,	"jpsikvf_ndf[n_jpsi]/F");  
  m_tree_jpsi->Branch("jpsi_L3D",	&m_jpsikvf_L3D,	"jpsikvf_L3D[n_jpsi]/F");  
  m_tree_jpsi->Branch("jpsi_sigmaL3D", &m_jpsikvf_sigmaL3D, "jpsikvf_sigmaL3D[n_jpsi]/F");  
  m_tree_jpsi->Branch("jpsi_L3DoverSigmaL3D", &m_jpsikvf_L3DoverSigmaL3D, "jpsikvf_L3DoverSigmaL3D[n_jpsi]/F");  

  m_tree_mujet = NULL;
  m_tree_mujet = new TTree("muTaggedJet_PF", "Mu-tagged jet info");  
  m_tree_mujet->Branch("n_mujet", &m_mujet_size, "n_mujet/I");
  m_tree_mujet->Branch("mujet_jet_btag_CSV", &m_mujet_jet_btag_CSV, "mujet_jet_btag_CSV[n_mujet]/F"); 
  m_tree_mujet->Branch("mujet_jet_4vector", "TClonesArray", &m_mujet_jet_lorentzvector, 1000, 0); 
  m_tree_mujet->Branch("mujet_jet_scaleFactor", &m_mujet_jet_scaleFactors.getBackingArray());
  m_tree_mujet->Branch("mujet_nonisomuplus_4vector", "TClonesArray", &m_mujet_nonisomuplus_lorentzvector, 1000, 0); 
  m_tree_mujet->Branch("mujet_nonisomuplus_pdgid", &m_mujet_nonisomuplus_pdgid, "mujet_nonisomuplus_pdgid[n_mujet]/I");
  for (auto& it: m_mujet_nonisomuplus_scaleFactors) {
    m_tree_mujet->Branch(("mujet_nonisomuplus_"+it.first).c_str(), & it.second.getBackingArray());
  }
  m_tree_mujet->Branch("mujet_nonisomuminus_4vector", "TClonesArray", &m_mujet_nonisomuminus_lorentzvector, 1000, 0); 
  m_tree_mujet->Branch("mujet_nonisomuminus_pdgid", &m_mujet_nonisomuminus_pdgid, "mujet_nonisomuminus_pdgid[n_mujet]/I");
  for (auto& it: m_mujet_nonisomuminus_scaleFactors) {
    m_tree_mujet->Branch(("mujet_nonisomuminus_"+it.first).c_str(), & it.second.getBackingArray());
  }
  m_tree_mujet->Branch("mujet_ntr", &m_mujet_ntr, "mujet_ntr[n_mujet]/I"); 
  m_tree_mujet->Branch("mujet_sump", &m_mujet_sump, "mujet_sump[n_mujet]/F"); 
  m_tree_mujet->Branch("mujet_sumpt", &m_mujet_sumpt, "mujet_sumpt[n_mujet]/F"); 
  m_tree_mujet->Branch("mujet_sumvecp", &m_mujet_sumvecp, "mujet_sumvecp[n_mujet]/F"); 
  m_tree_mujet->Branch("n_tr", &m_mujet_tr_size, "n_tr/I");
  m_tree_mujet->Branch("mujet_tr_indmujet", &m_mujet_tr_indmujet, "mujet_tr_indmujet[n_tr]/I"); 
  m_tree_mujet->Branch("mujet_tr_pdgid", &m_mujet_tr_pdgid, "mujet_tr_pdgid[n_tr]/I"); 
  m_tree_mujet->Branch("mujet_tr_4vector", "TClonesArray", &m_mujet_tr_lorentzvector, 1000, 0); 
  m_tree_mujet->Branch("n_d0", &m_mujet_d0_size, "n_d0/I");
  m_tree_mujet->Branch("mujet_nd0", &m_mujet_nd0, "mujet_nd0[n_mujet]/I");
  m_tree_mujet->Branch("mujet_d0_indmujet", &m_mujet_d0kvf_indmujet, "mujet_d0_indmujet[n_d0]/I"); 
  m_tree_mujet->Branch("mujet_d0pf_4vector", "TClonesArray", &m_mujet_d0pf_lorentzvector, 1000, 0);
  m_tree_mujet->Branch("mujet_d0_4vector", "TClonesArray", &m_mujet_d0kvf_lorentzvector, 1000, 0);
  m_tree_mujet->Branch("mujet_d0_kaon_4vector", "TClonesArray", &m_mujet_d0kvf_kaon_lorentzvector, 1000, 0);
  m_tree_mujet->Branch("mujet_d0_kaon_pdgid", &m_mujet_d0kvf_kaon_pdgid, "mujet_d0_kaon_pdgid[n_d0]/I"); 
  m_tree_mujet->Branch("mujet_d0_pion_4vector", "TClonesArray", &m_mujet_d0kvf_pion_lorentzvector, 1000, 0);
  m_tree_mujet->Branch("mujet_d0_pion_pdgid", &m_mujet_d0kvf_pion_pdgid, "mujet_d0_pion_pdgid[n_d0]/I"); 
  m_tree_mujet->Branch("mujet_d0_vx", &m_mujet_d0kvf_vx,	"mujet_d0kvf_vx[n_d0]/F");  
  m_tree_mujet->Branch("mujet_d0_vy", &m_mujet_d0kvf_vy,	"mujet_d0kvf_vy[n_d0]/F");  
  m_tree_mujet->Branch("mujet_d0_vz", &m_mujet_d0kvf_vz,	"mujet_d0kvf_vz[n_d0]/F");
  m_tree_mujet->Branch("mujet_d0_vtxvalid", &m_mujet_d0kvf_vtxvalid, "mujet_d0kvf_vtxvalid[n_d0]/O");
  m_tree_mujet->Branch("mujet_d0_vtxchi2", &m_mujet_d0kvf_vtxchi2,  "mujet_d0kvf_vtxchi2[n_d0]/F");
  m_tree_mujet->Branch("mujet_d0_ndf",	&m_mujet_d0kvf_ndf,	"mujet_d0kvf_ndf[n_d0]/F");  
  m_tree_mujet->Branch("mujet_d0_L3D",	&m_mujet_d0kvf_L3D,	"mujet_d0kvf_L3D[n_d0]/F");  
  m_tree_mujet->Branch("mujet_d0_sigmaL3D", &m_mujet_d0kvf_sigmaL3D, "mujet_d0kvf_sigmaL3D[n_d0]/F");  
  m_tree_mujet->Branch("mujet_d0_L3DoverSigmaL3D", &m_mujet_d0kvf_L3DoverSigmaL3D, "mujet_d0kvf_L3DoverSigmaL3D[n_d0]/F");  
  m_tree_mujet->Branch("n_unfold_tr", &m_mujet_unfold_tr_size, "n_unfold_tr/I");
  m_tree_mujet->Branch("mujet_unfold_indmujet", &m_mujet_unfold_indmujet, "mujet_unfold_indmujet[n_unfold_tr]/I"); 
  m_tree_mujet->Branch("mujet_unfold_tr_recopt", &m_mujet_unfold_tr_recopt, "mujet_unfold_tr_recopt[n_unfold_tr]/F");
  m_tree_mujet->Branch("mujet_unfold_tr_recoeta", &m_mujet_unfold_tr_recoeta, "mujet_unfold_tr_recoeta[n_unfold_tr]/F");
  m_tree_mujet->Branch("mujet_unfold_tr_genpt", &m_mujet_unfold_tr_genpt, "mujet_unfold_tr_genpt[n_unfold_tr]/F");
  m_tree_mujet->Branch("mujet_unfold_tr_geneta", &m_mujet_unfold_tr_geneta, "mujet_unfold_tr_geneta[n_unfold_tr]/F");
  m_tree_mujet->Branch("mujet_unfold_tr_dr", &m_mujet_unfold_tr_dr, "mujet_unfold_tr_dr[n_unfold_tr]/F");
  m_tree_mujet->Branch("mujet_unfold_mu_recopt", &m_mujet_unfold_mu_recopt, "mujet_unfold_mu_recopt[n_unfold_tr]/F");
  m_tree_mujet->Branch("mujet_unfold_mu_recoeta", &m_mujet_unfold_mu_recoeta, "mujet_unfold_mu_recoeta[n_unfold_tr]/F");
  m_tree_mujet->Branch("mujet_unfold_mu_genpt", &m_mujet_unfold_mu_genpt, "mujet_unfold_mu_genpt[n_unfold_tr]/F");
  m_tree_mujet->Branch("mujet_unfold_mu_geneta", &m_mujet_unfold_mu_geneta, "mujet_unfold_mu_geneta[n_unfold_tr]/F");
  m_tree_mujet->Branch("mujet_unfold_mu_dr", &m_mujet_unfold_mu_dr, "mujet_unfold_mu_dr[n_unfold_tr]/F");
}

void KVFExtractor::beginJob() {
  if (!mUseGlobalTagForJEC) {
    mTxtCorrector = makeFactorizedJetCorrectorFromXML(mJecPayload, mJecJetAlgo, m_isMC);
    std::cout << "Using text files for JEC" << std::endl;
  } else {
    std::cout << "Using global tag for JEC" << std::endl;
  }
}


  KVFExtractor::KVFExtractor(const std::string& name_jpsi, const std::string& name_d0, const std::string& name_jet, std::shared_ptr<ScaleFactorService> sf, TFile *a_file)
: BaseExtractor(name_jet, sf)
{

  std::cout << "KVFExtractor objet is retrieved" << std::endl;
  m_file = a_file;

  // Tree definition
  m_OK = false;

  m_tree_jpsi = dynamic_cast<TTree*>(a_file->Get(name_jpsi.c_str()));
  m_tree_mujet = dynamic_cast<TTree*>(a_file->Get("muTaggedJet_PF"));
  const auto& sfWorkingPoints = m_scaleFactorService->getMuonScaleFactorWorkingPoints();

  if (m_tree_jpsi) {

    m_jpsi_jet_lorentzvector    = new TClonesArray("TLorentzVector");
    m_jpsi_jet_scaleFactors.setWriteMode();
    m_jpsipf_lorentzvector      = new TClonesArray("TLorentzVector");
    m_jpsikvf_lorentzvector     = new TClonesArray("TLorentzVector");
    m_jpsikvf_mu1_lorentzvector = new TClonesArray("TLorentzVector");
    for (auto& it: sfWorkingPoints) {
      std::string name = "muon_scaleFactor_" + ScaleFactorService::workingPointToString(it.first) + "eff_" + ScaleFactorService::workingPointToString(it.second) + "iso";
      m_jpsikvf_mu1_muon_scaleFactors[name] = ScaleFactorCollection();
      m_jpsikvf_mu1_muon_scaleFactors[name].setWriteMode();
    }
    m_jpsikvf_mu2_lorentzvector = new TClonesArray("TLorentzVector");
    for (auto& it: sfWorkingPoints) {
      std::string name = "muon_scaleFactor_" + ScaleFactorService::workingPointToString(it.first) + "eff_" + ScaleFactorService::workingPointToString(it.second) + "iso";
      m_jpsikvf_mu2_muon_scaleFactors[name] = ScaleFactorCollection();
      m_jpsikvf_mu2_muon_scaleFactors[name].setWriteMode();
    }

    if (m_tree_jpsi->FindBranch("n_jpsi")) 
      m_tree_jpsi->Branch("n_jpsi", &m_jpsi_size);
    if (m_tree_jpsi->FindBranch("jpsi_indjet")) 
      m_tree_jpsi->Branch("jpsi_indjet", &m_jpsi_indjet);
    if (m_tree_jpsi->FindBranch("jpsi_jet_btag_CSV")) 
      m_tree_jpsi->Branch("jpsi_jet_btag_CSV", &m_jpsi_jet_btag_CSV);
    if (m_tree_jpsi->FindBranch("jpsi_jet_4vector")) 
      m_tree_jpsi->Branch("jpsi_jet_4vector", &m_jpsi_jet_lorentzvector);
    if (m_tree_jpsi->FindBranch("jpsi_jet_scaleFactor"))
      m_tree_jpsi->SetBranchAddress("jpsi_jet_scaleFactor", &m_jpsi_jet_scaleFactors.getBackingArray());
    if (m_tree_jpsi->FindBranch("jpsi_indpf1")) 
      m_tree_jpsi->Branch("jpsi_indpf1", &m_jpsi_indpf1);
    if (m_tree_jpsi->FindBranch("jpsi_indpf2")) 
      m_tree_jpsi->Branch("jpsi_indpf2", &m_jpsi_indpf2);
    if (m_tree_jpsi->FindBranch("jpsipf_4vector")) 
      m_tree_jpsi->Branch("jpsipf_4vector", &m_jpsipf_lorentzvector);
    if (m_tree_jpsi->FindBranch("jpsi_4vector")) 
      m_tree_jpsi->Branch("jpsi_4vector", &m_jpsikvf_lorentzvector);
    if (m_tree_jpsi->FindBranch("jpsi_mu1_4vector")) 
      m_tree_jpsi->Branch("jpsi_mu1_4vector", &m_jpsikvf_mu1_lorentzvector);
    for (auto& it: m_jpsikvf_mu1_muon_scaleFactors) {
      if (m_tree_jpsi->FindBranch(("jpsi_mu1_"+it.first).c_str()))
        m_tree_jpsi->Branch(("jpsi_mu1_"+it.first).c_str(), & it.second.getBackingArray());
    }
    if (m_tree_jpsi->FindBranch("jpsi_mu2_4vector")) 
      m_tree_jpsi->Branch("jpsi_mu2_4vector", &m_jpsikvf_mu2_lorentzvector);
    for (auto& it: m_jpsikvf_mu2_muon_scaleFactors) {
      if (m_tree_jpsi->FindBranch(("jpsi_mu2_"+it.first).c_str()))
        m_tree_jpsi->Branch(("jpsi_mu2_"+it.first).c_str(), & it.second.getBackingArray());
    }
    if (m_tree_jpsi->FindBranch("jpsi_vx")) 
      m_tree_jpsi->Branch("jpsi_vx", &m_jpsikvf_vx);  
    if (m_tree_jpsi->FindBranch("jpsi_vy")) 
      m_tree_jpsi->Branch("jpsi_vy", &m_jpsikvf_vy);  
    if (m_tree_jpsi->FindBranch("jpsi_vz")) 
      m_tree_jpsi->Branch("jpsi_vz", &m_jpsikvf_vz);
    if (m_tree_jpsi->FindBranch("jpsi_vtxvalid")) 
      m_tree_jpsi->Branch("jpsi_vtxvalid", &m_jpsikvf_vtxvalid);
    if (m_tree_jpsi->FindBranch("jpsi_vtxchi2")) 
      m_tree_jpsi->Branch("jpsi_vtxchi2", &m_jpsikvf_vtxchi2);
    if (m_tree_jpsi->FindBranch("jpsi_ndf")) 
      m_tree_jpsi->Branch("jpsi_ndf", &m_jpsikvf_ndf);  
    if (m_tree_jpsi->FindBranch("jpsi_L3D")) 
      m_tree_jpsi->Branch("jpsi_L3D", &m_jpsikvf_L3D);  
    if (m_tree_jpsi->FindBranch("jpsi_sigmaL3D")) 
      m_tree_jpsi->Branch("jpsi_sigmaL3D", &m_jpsikvf_sigmaL3D);  
    if (m_tree_jpsi->FindBranch("jpsi_L3DoverSigmaL3D")) 
      m_tree_jpsi->Branch("jpsi_L3DoverSigmaL3D", &m_jpsikvf_L3DoverSigmaL3D);  
  }
  if (m_tree_mujet) {

    m_mujet_jet_lorentzvector           = new TClonesArray("TLorentzVector");
    m_mujet_jet_scaleFactors.setWriteMode();
    m_mujet_nonisomuplus_lorentzvector  = new TClonesArray("TLorentzVector");
    for (auto& it: sfWorkingPoints) {
      std::string name = "muon_scaleFactor_" + ScaleFactorService::workingPointToString(it.first) + "eff_" + ScaleFactorService::workingPointToString(it.second) + "iso";
      m_mujet_nonisomuplus_scaleFactors[name] = ScaleFactorCollection();
      m_mujet_nonisomuplus_scaleFactors[name].setWriteMode();
    }
    m_mujet_nonisomuminus_lorentzvector = new TClonesArray("TLorentzVector");
    for (auto& it: sfWorkingPoints) {
      std::string name = "muon_scaleFactor_" + ScaleFactorService::workingPointToString(it.first) + "eff_" + ScaleFactorService::workingPointToString(it.second) + "iso";
      m_mujet_nonisomuminus_scaleFactors[name] = ScaleFactorCollection();
      m_mujet_nonisomuminus_scaleFactors[name].setWriteMode();
    }
    m_mujet_tr_lorentzvector            = new TClonesArray("TLorentzVector");
    m_mujet_d0pf_lorentzvector          = new TClonesArray("TLorentzVector");
    m_mujet_d0kvf_lorentzvector         = new TClonesArray("TLorentzVector");
    m_mujet_d0kvf_pion_lorentzvector    = new TClonesArray("TLorentzVector");
    m_mujet_d0kvf_kaon_lorentzvector    = new TClonesArray("TLorentzVector");

    if (m_tree_mujet->FindBranch("n_mujet")) 
      m_tree_mujet->Branch("n_mujet", &m_mujet_size);
    if (m_tree_mujet->FindBranch("mujet_jet_btag_CSV")) 
      m_tree_mujet->Branch("mujet_jet_btag_CSV", &m_mujet_jet_btag_CSV);
    if (m_tree_mujet->FindBranch("mujet_jet_4vector")) 
      m_tree_mujet->Branch("mujet_jet_4vector", &m_mujet_jet_lorentzvector);
    if (m_tree_mujet->FindBranch("mujet_jet_scaleFactor"))
      m_tree_mujet->SetBranchAddress("mujet_jet_scaleFactor", &m_mujet_jet_scaleFactors.getBackingArray());
    if (m_tree_mujet->FindBranch("mujet_nonisomuplus_4vector")) 
      m_tree_mujet->Branch("mujet_nonisomuplus_4vector", &m_mujet_nonisomuplus_lorentzvector);
    if (m_tree_mujet->FindBranch("mujet_nonisomuplus_pdgid")) 
      m_tree_mujet->Branch("mujet_nonisomuplus_pdgid", &m_mujet_nonisomuplus_pdgid);
    for (auto& it: m_mujet_nonisomuplus_scaleFactors) {
      if (m_tree_mujet->FindBranch(("mujet_nonisomuplus_"+it.first).c_str()))
        m_tree_mujet->Branch(("mujet_nonisomuplus_"+it.first).c_str(), & it.second.getBackingArray());
    }
    if (m_tree_mujet->FindBranch("mujet_nonisomuminus_4vector")) 
      m_tree_mujet->Branch("mujet_nonisomuminus_4vector", &m_mujet_nonisomuminus_lorentzvector);
    if (m_tree_mujet->FindBranch("mujet_nonisomuminus_pdgid")) 
      m_tree_mujet->Branch("mujet_nonisomuminus_pdgid", &m_mujet_nonisomuminus_pdgid);
    for (auto& it: m_mujet_nonisomuminus_scaleFactors) {
      if (m_tree_mujet->FindBranch(("mujet_nonisomuminus_"+it.first).c_str()))
        m_tree_mujet->Branch(("mujet_nonisomuminus_"+it.first).c_str(), & it.second.getBackingArray());
    }
    if (m_tree_mujet->FindBranch("mujet_ntr")) 
      m_tree_mujet->Branch("mujet_ntr", &m_mujet_ntr);
    if (m_tree_mujet->FindBranch("mujet_sump")) 
      m_tree_mujet->Branch("mujet_sump", &m_mujet_sump);
    if (m_tree_mujet->FindBranch("mujet_sumpt")) 
      m_tree_mujet->Branch("mujet_sumpt", &m_mujet_sumpt);
    if (m_tree_mujet->FindBranch("mujet_sumvecp")) 
      m_tree_mujet->Branch("mujet_sumvecp", &m_mujet_sumvecp);
    if (m_tree_mujet->FindBranch("n_tr")) 
      m_tree_mujet->Branch("n_tr", &m_mujet_tr_size);
    if (m_tree_mujet->FindBranch("mujet_tr_indmujet")) 
      m_tree_mujet->Branch("mujet_tr_indmujet", &m_mujet_tr_indmujet); 
    if (m_tree_mujet->FindBranch("mujet_tr_pdgid")) 
      m_tree_mujet->Branch("mujet_tr_pdgid", &m_mujet_tr_pdgid); 
    if (m_tree_mujet->FindBranch("mujet_tr_4vector")) 
      m_tree_mujet->Branch("mujet_tr_4vector", &m_mujet_tr_lorentzvector); 
    if (m_tree_mujet->FindBranch("n_d0")) 
      m_tree_mujet->Branch("n_d0", &m_mujet_d0_size);
    if (m_tree_mujet->FindBranch("mujet_nd0")) 
      m_tree_mujet->Branch("mujet_nd0", &m_mujet_nd0);
    if (m_tree_mujet->FindBranch("mujet_d0_indmujet")) 
      m_tree_mujet->Branch("mujet_d0_indmujet", &m_mujet_d0kvf_indmujet); 
    if (m_tree_mujet->FindBranch("mujet_d0pf_4vector")) 
      m_tree_mujet->Branch("mujet_d0pf_4vector", &m_mujet_d0pf_lorentzvector);
    if (m_tree_mujet->FindBranch("mujet_d0_4vector")) 
      m_tree_mujet->Branch("mujet_d0_4vector", &m_mujet_d0kvf_lorentzvector);
    if (m_tree_mujet->FindBranch("mujet_d0_kaon_4vector")) 
      m_tree_mujet->Branch("mujet_d0_kaon_4vector", &m_mujet_d0kvf_kaon_lorentzvector);
    if (m_tree_mujet->FindBranch("mujet_d0_kaon_pdgid")) 
      m_tree_mujet->Branch("mujet_d0_kaon_pdgid", &m_mujet_d0kvf_kaon_pdgid); 
    if (m_tree_mujet->FindBranch("mujet_d0_pion_4vector")) 
      m_tree_mujet->Branch("mujet_d0_pion_4vector", &m_mujet_d0kvf_pion_lorentzvector);
    if (m_tree_mujet->FindBranch("mujet_d0_pion_pdgid")) 
      m_tree_mujet->Branch("mujet_d0_pion_pdgid", &m_mujet_d0kvf_pion_pdgid); 
    if (m_tree_mujet->FindBranch("mujet_d0_vx")) 
      m_tree_mujet->Branch("mujet_d0_vx", &m_mujet_d0kvf_vx);  
    if (m_tree_mujet->FindBranch("mujet_d0_vy")) 
      m_tree_mujet->Branch("mujet_d0_vy", &m_mujet_d0kvf_vy);  
    if (m_tree_mujet->FindBranch("mujet_d0_vz")) 
      m_tree_mujet->Branch("mujet_d0_vz", &m_mujet_d0kvf_vz);
    if (m_tree_mujet->FindBranch("mujet_d0_vtxvalid")) 
      m_tree_mujet->Branch("mujet_d0_vtxvalid", &m_mujet_d0kvf_vtxvalid);
    if (m_tree_mujet->FindBranch("mujet_d0_vtxchi2")) 
      m_tree_mujet->Branch("mujet_d0_vtxchi2", &m_mujet_d0kvf_vtxchi2);
    if (m_tree_mujet->FindBranch("mujet_d0_ndf")) 
      m_tree_mujet->Branch("mujet_d0_ndf",	&m_mujet_d0kvf_ndf);  
    if (m_tree_mujet->FindBranch("mujet_d0_L3D")) 
      m_tree_mujet->Branch("mujet_d0_L3D",	&m_mujet_d0kvf_L3D);  
    if (m_tree_mujet->FindBranch("mujet_d0_sigmaL3D")) 
      m_tree_mujet->Branch("mujet_d0_sigmaL3D", &m_mujet_d0kvf_sigmaL3D);  
    if (m_tree_mujet->FindBranch("mujet_d0_L3DoverSigmaL3D")) 
      m_tree_mujet->Branch("mujet_d0_L3DoverSigmaL3D", &m_mujet_d0kvf_L3DoverSigmaL3D);  
    if (m_tree_mujet->FindBranch("n_unfold_tr")) 
      m_tree_mujet->Branch("n_unfold_tr", &m_mujet_unfold_tr_size);
    if (m_tree_mujet->FindBranch("mujet_unfold_indmujet")) 
      m_tree_mujet->Branch("mujet_unfold_indmujet", &m_mujet_unfold_indmujet); 
    if (m_tree_mujet->FindBranch("mujet_unfold_tr_recopt")) 
      m_tree_mujet->Branch("mujet_unfold_tr_recopt", &m_mujet_unfold_tr_recopt);
    if (m_tree_mujet->FindBranch("mujet_unfold_tr_recoeta")) 
      m_tree_mujet->Branch("mujet_unfold_tr_recoeta", &m_mujet_unfold_tr_recoeta);
    if (m_tree_mujet->FindBranch("mujet_unfold_tr_genpt")) 
      m_tree_mujet->Branch("mujet_unfold_tr_genpt", &m_mujet_unfold_tr_genpt);
    if (m_tree_mujet->FindBranch("mujet_unfold_tr_geneta")) 
      m_tree_mujet->Branch("mujet_unfold_tr_geneta", &m_mujet_unfold_tr_geneta);
    if (m_tree_mujet->FindBranch("mujet_unfold_tr_dr")) 
      m_tree_mujet->Branch("mujet_unfold_tr_dr", &m_mujet_unfold_tr_dr);
    if (m_tree_mujet->FindBranch("mujet_unfold_mu_recopt")) 
      m_tree_mujet->Branch("mujet_unfold_mu_recopt", &m_mujet_unfold_mu_recopt);
    if (m_tree_mujet->FindBranch("mujet_unfold_mu_recoeta")) 
      m_tree_mujet->Branch("mujet_unfold_mu_recoeta", &m_mujet_unfold_mu_recoeta);
    if (m_tree_mujet->FindBranch("mujet_unfold_mu_genpt")) 
      m_tree_mujet->Branch("mujet_unfold_mu_genpt", &m_mujet_unfold_mu_genpt);
    if (m_tree_mujet->FindBranch("mujet_unfold_mu_geneta")) 
      m_tree_mujet->Branch("mujet_unfold_mu_geneta", &m_mujet_unfold_mu_geneta);
    if (m_tree_mujet->FindBranch("mujet_unfold_mu_dr")) 
      m_tree_mujet->Branch("mujet_unfold_mu_dr", &m_mujet_unfold_mu_dr);
  }

  m_OK = true;

}


KVFExtractor::~KVFExtractor()
{
  delete mTxtCorrector;
}


bool KVFExtractor::isPFJetLoose(const pat::Jet& jet)
{
  // Taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetID

  if (! jet.isPFJet())
    return false;

  // Jet ID works only with uncorrected jet. *EnergyFraction functions take care of that all alone
  bool isValid = jet.neutralHadronEnergyFraction() < 0.99;
  isValid &= jet.neutralEmEnergyFraction() < 0.99;
  isValid &= jet.getPFConstituents().size() > 1;
  if (fabs(jet.eta()) < 2.4) {
    isValid &= jet.chargedHadronEnergyFraction() > 0.;
    isValid &= jet.chargedMultiplicity() > 0;
    isValid &= jet.chargedEmEnergyFraction() < 0.99;
  }

  return isValid;
}


//
// Method filling the main particle tree
//

void KVFExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* m_MC) 
{
  edm::Handle<pat::JetCollection>  jetHandle;
  event.getByLabel(m_tag, jetHandle);
  pat::JetCollection p_jets = *jetHandle;

  edm::Handle<std::vector<reco::Vertex>>  pvHandle;
  event.getByLabel(m_vertexTag, pvHandle);

  bool unfold = false;
  edm::Handle<reco::GenParticleCollection> genParticles;
  if (event.getByLabel("genParticles", genParticles)) unfold = true;

  reset();
  m_size = 0;

  if (mJESSign != 0) {
    if (! jecUncertainty.get()) {
      if (mJecFilename.length() > 0) {
        std::cout << "Reading JES uncertainties from '" << mJecFilename << "'" << std::endl;
        jecUncertainty.reset(new JetCorrectionUncertainty(mJecFilename));
      } else {
        std::cout << "Reading JES uncertainties from Global Tag" << std::endl;
        edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
        iSetup.get<JetCorrectionsRecord>().get("AK5PFchs",JetCorParColl);
        JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
        jecUncertainty.reset(new JetCorrectionUncertainty(JetCorPar));
      }
    }
  }

  if (mCorrectJets) {
    correctJets(p_jets, event, iSetup);
  } else {
    extractRawJets(p_jets);
  }

  //do Jets MET resolution corrections
  if (m_MC && mDoJER)
    correctJetsResolution(p_jets);

  // JES systematics
  if (mJESSign != 0)
    doJESSystematics(p_jets);

  int nJpsi = 0; 
  int nMuJet = 0;
  int nTr = 0;
  int nUnfoldTr = 0;
  int nD0 = 0;

  for (unsigned int i = 0; i < p_jets.size(); ++i)
  {
    const pat::Jet& rawJet = *(p_jets.at(i).userData<pat::Jet>("rawJet"));
#if DEBUG
    if (! mCorrectJets) {
      // Test if jet id works the same on our stored raw jets and the jet itself
      if (isPFJetLoose(rawJet) != isPFJetLoose(p_jets.at(i))) {
        std::cout << "Error: there's something wrong with your jets!" << std::endl;
      }
    }
#endif
    if(mDoLooseJetID) {
      if (! isPFJetLoose(rawJet)) 
        continue;
    }


    pat::JetRef jetRef(jetHandle, i); 
    bool hasNonIsoMu = false;
    float SumP = 0.;
    float SumPt = 0.;
    TLorentzVector SumVecP;
    SumVecP.SetPxPyPzE(0., 0., 0., 0.); 
    std::vector<reco::PFCandidate> myPFs; 
    std::vector<reco::PFCandidate> myKPis; 
    std::vector<reco::PFCandidate> myPFs2Unfold; 


    // Reconstruct the J/psi
    //----------------------

    const std::vector<reco::PFCandidatePtr>& PFpart = p_jets.at(i).getPFConstituents();
    unsigned int npfs = PFpart.size(); 

    ParticleMass muon_mass  = 0.1056583;
    float        muon_sigma = 0.0000001;

    // ParticleMass jpsi_mass  = 3.09687;
    // float        jpsi_sigma = 0.00009;

    // To transform Track to TransientTrack, first need to get the builder:
    edm::ESHandle<TransientTrackBuilder> theMuB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theMuB);  

    for(unsigned int j = 0; j < npfs; ++j) {
      if (!PFpart[j]->trackRef()) continue;

      if (PFpart[j]->trackRef()->quality(reco::Track::tight)) {
        if (PFpart[j]->pt() > m_trUnfoldMinPt)
          myPFs2Unfold.push_back(*PFpart[j]);
        if (PFpart[j]->pt() > m_trSumMinPt) {
          SumP += PFpart[j]->p();
          SumPt += PFpart[j]->pt();
          TLorentzVector VecP;
          VecP.SetPxPyPzE(PFpart[j]->px(), PFpart[j]->py(), PFpart[j]->pz(), PFpart[j]->energy());
          SumVecP = SumVecP + VecP;  
          myPFs.push_back(*PFpart[j]);
          if (abs(PFpart[j]->pdgId()) == 13) hasNonIsoMu = true;
          else myKPis.push_back(*PFpart[j]);
        }
      } // compute some useful variables for mu tagged jets

      // Both PF particle should be a muon
      /*
      if (abs(PFpart[j]->pdgId()) != 13) continue;
      if (PFpart[j]->pt() < m_muJpsiMinPt) continue;
      */
      // try BPH selection https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Soft_Muon
      if (!PFpart[j]->muonRef()) continue;
      if (!muon::isGoodMuon(*(PFpart[j]->muonRef()), muon::TMOneStationTight)) continue;
      if (PFpart[j]->muonRef()->innerTrack()->hitPattern().trackerLayersWithMeasurement() <= 5) continue;
      if (PFpart[j]->muonRef()->innerTrack()->hitPattern().pixelLayersWithMeasurement() <= 0) continue;
      if (!PFpart[j]->muonRef()->innerTrack()->quality(reco::TrackBase::highPurity)) continue;
      if (fabs(PFpart[j]->muonRef()->innerTrack()->dxy(pvHandle->at(0).position())) >= 0.3 || fabs(PFpart[j]->muonRef()->innerTrack()->dz(pvHandle->at(0).position())) >= 20.) continue;

      for(unsigned int k = j+1; k < npfs; ++k) {

        // Both PF particle should be a muon
        /*
        if (abs(PFpart[k]->pdgId()) != 13) continue;
        if (PFpart[k]->pt() < m_muJpsiMinPt) continue;
        */
        // try BPH selection https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Soft_Muon
        if (!PFpart[k]->muonRef()) continue;
        if (!muon::isGoodMuon(*(PFpart[k]->muonRef()), muon::TMOneStationTight)) continue;
        if (PFpart[k]->muonRef()->innerTrack()->hitPattern().trackerLayersWithMeasurement() <= 5) continue;
        if (PFpart[k]->muonRef()->innerTrack()->hitPattern().pixelLayersWithMeasurement() <= 0) continue;
        if (!PFpart[k]->muonRef()->innerTrack()->quality(reco::TrackBase::highPurity)) continue;
        if (fabs(PFpart[k]->muonRef()->innerTrack()->dxy(pvHandle->at(0).position())) >= 0.3 || fabs(PFpart[k]->muonRef()->innerTrack()->dz(pvHandle->at(0).position())) >= 20.) continue;

        // Both PF particle should be of OS
        if (!PFpart[k]->trackRef()) continue;
        if (PFpart[j]->charge() + PFpart[k]->charge() != 0) continue;

        double eJpsi  = PFpart[j]->energy()+PFpart[k]->energy();
        double pxJpsi = PFpart[j]->px()+PFpart[k]->px();
        double pyJpsi = PFpart[j]->py()+PFpart[k]->py();
        double pzJpsi = PFpart[j]->pz()+PFpart[k]->pz();
        double mJpsi = pow(eJpsi,2)-pow(pxJpsi,2)-pow(pyJpsi,2)-pow(pzJpsi,2);
        if (mJpsi > 0.) mJpsi = sqrt(mJpsi);
        else mJpsi = 0.;

        if (mJpsi >= m_jpsiMassMin && mJpsi < m_jpsiMassMax) {
          if (nJpsi > 0) std::cout << "Warning: another Jpsi found with this PF particle" << std::endl;
          ++nJpsi;

          // Fill tree info for links between Jpsi and jet
          m_jpsi_indjet[nJpsi-1] = i;
          m_jpsi_jet_btag_CSV[nJpsi-1] = p_jets.at(i).bDiscriminator("combinedSecondaryVertexBJetTags");
          new((*m_jpsi_jet_lorentzvector)[nJpsi-1]) TLorentzVector((p_jets.at(i)).px(),(p_jets.at(i)).py(),(p_jets.at(i)).pz(),(p_jets.at(i)).energy());
          if (m_isMC) {
            int mcFlavor = abs(p_jets.at(i).partonFlavour());
            ScaleFactorService::Flavor flavor = ScaleFactorService::B;
            if (mcFlavor == 4) {
              flavor = ScaleFactorService::C;
            } else if ((mcFlavor <= 3) || (mcFlavor == 21)) {
              // If mcFlavor == 0, assume it's a light jet
              flavor = ScaleFactorService::LIGHT;
            }
            m_jpsi_jet_scaleFactors.push_back(m_scaleFactorService->getBTaggingScaleFactor(flavor, p_jets.at(i).et(), p_jets.at(i).eta()));
          }
          // Fill tree info for links between Jpsi and PF particles
          m_jpsi_indpf1[nJpsi-1] = j;
          m_jpsi_indpf2[nJpsi-1] = k;
          // save the Jpsi from 2 PF muons (which can be different from Jpsi from 2 tracks)
          new((*m_jpsipf_lorentzvector)[nJpsi-1]) TLorentzVector(pxJpsi, pyJpsi, pzJpsi, eJpsi);

          // Make the Transient tracks
          reco::TransientTrack tr1 = (*theMuB).build(PFpart[j]->trackRef());
          reco::TransientTrack tr2 = (*theMuB).build(PFpart[k]->trackRef());

          // A complicated Kalman fit :
          //---------------------------

          //Creating a KinematicParticleFactory
          KinematicParticleFactoryFromTransientTrack pFactory;

          //initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered
          float chi = 0.;
          float ndf = 0.;

          //making particles
          std::vector<RefCountedKinematicParticle> muonParticles;
          muonParticles.push_back(pFactory.particle (tr1,muon_mass,chi,ndf,muon_sigma));
          muonParticles.push_back(pFactory.particle (tr2,muon_mass,chi,ndf,muon_sigma));

          /* Example of a simple vertex fit, without other constraints
           * The reconstructed decay tree is a result of the kinematic fit
           * The KinematicParticleVertexFitter fits the final state particles to their vertex and
           * reconstructs the decayed state
           */

          // creating the vertex fitter
          KinematicParticleVertexFitter fitter;
          // reconstructing a J/Psi decay
          RefCountedKinematicTree vertexFitTree = fitter.fit(muonParticles);

          if (!vertexFitTree->isValid()) {
            std::cout <<"J/psi vertexTree is invalid. Fit failed." << std::endl;

            // Need to fill empty quantities : OK for tables, but need to create empty TLorentzVector
            new((*m_jpsikvf_lorentzvector)[nJpsi-1]) TLorentzVector(0,0,0,0);
            new((*m_jpsikvf_mu1_lorentzvector)[nJpsi-1]) TLorentzVector(0,0,0,0);
            new((*m_jpsikvf_mu2_lorentzvector)[nJpsi-1]) TLorentzVector(0,0,0,0);

            continue;

          } else {
            //accessing the tree components, move pointer to top
            vertexFitTree->movePointerToTheTop();

            //We are now at the top of the decay tree getting the jpsi reconstructed KinematicPartlcle
            RefCountedKinematicParticle jpsi1 = vertexFitTree->currentParticle();
            AlgebraicVector7 par0 = jpsi1->currentState().kinematicParameters().vector();
            double e0 = jpsi1->currentState().kinematicParameters().energy();
            new((*m_jpsikvf_lorentzvector)[nJpsi-1]) TLorentzVector(par0(3),par0(4),par0(5),e0);

            RefCountedKinematicVertex jpsi1_vertex = vertexFitTree->currentDecayVertex();
            if ( jpsi1_vertex->vertexIsValid()) {

              m_jpsikvf_vx[nJpsi-1] = jpsi1_vertex->position().x();
              m_jpsikvf_vy[nJpsi-1] = jpsi1_vertex->position().y();
              m_jpsikvf_vz[nJpsi-1] = jpsi1_vertex->position().z();
              m_jpsikvf_vtxvalid[nJpsi-1] = true;
              m_jpsikvf_vtxchi2[nJpsi-1] = jpsi1_vertex->chiSquared();
              m_jpsikvf_ndf[nJpsi-1] = jpsi1_vertex->degreesOfFreedom();

              // Compute the distance between the PV and the Jpsi vertex :
              //----------------------------------------------------------

              const reco::VertexCollection vtx = *(pvHandle.product());

              GlobalPoint svPos    = jpsi1_vertex->position();
              GlobalError svPosErr = jpsi1_vertex->error();

              // If the  PV does not exist, compute distance wrt to the detector center (0,0,0)

              double sigmax = 0.;
              double sigmay = 0.;
              double sigmaz = 0.;
              if (vtx.size() > 0) {
                sigmax = sqrt(vtx[0].xError()*vtx[0].xError() + svPosErr.cxx()*svPosErr.cxx());
                sigmay = sqrt(vtx[0].yError()*vtx[0].yError() + svPosErr.cyy()*svPosErr.cyy());
                sigmaz = sqrt(vtx[0].zError()*vtx[0].zError() + svPosErr.czz()*svPosErr.czz());
              }  else {
                sigmax = sqrt(svPosErr.cxx()*svPosErr.cxx());
                sigmay = sqrt(svPosErr.cyy()*svPosErr.cyy());
                sigmaz = sqrt(svPosErr.czz()*svPosErr.czz());
              }

              double px  = par0(3);
              double py  = par0(4);
              double pz  = par0(5);
              double nrj = e0;
              double m = sqrt( nrj*nrj - px*px - py*py - pz*pz );

              double interx = pow((px/m)/sigmax, 2.);
              double intery = pow((py/m)/sigmay, 2.);
              double interz = pow((pz/m)/sigmaz, 2.);

              m_jpsikvf_sigmaL3D[nJpsi-1] = pow( interx + intery + interz , -0.5);
              //std::cout << "sigmaL3D = " << m_jpsikvf_sigmaL3D[nJpsi-1] << std::endl;

              double part1 = 0.;
              double part2 = 0.;
              double part3 = 0.;
              if (vtx.size() >0) {
                part1 = (px/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmax,2.)*(svPos.x() - vtx[0].x());
                part2 = (py/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmay,2.)*(svPos.y() - vtx[0].y());
                part3 = (pz/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmaz,2.)*(svPos.z() - vtx[0].z());
              }  else {
                part1 = (px/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmax,2.)*svPos.x();
                part2 = (py/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmay,2.)*svPos.y();
                part3 = (pz/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmaz,2.)*svPos.z();
              }

              m_jpsikvf_L3D[nJpsi-1] = fabs(part1 + part2 + part3);
              //std::cout << "L3D = " << m_jpsikvf_L3D[nJpsi-1] << std::endl;  

              m_jpsikvf_L3DoverSigmaL3D[nJpsi-1] = m_jpsikvf_L3D[nJpsi-1]/m_jpsikvf_sigmaL3D[nJpsi-1];
              //std::cout << "(L/sigma)3D = " << m_jpsikvf_L3DoverSigmaL3D[nJpsi-1] << std::endl;    

              if (!pvHandle.isValid()) {
                std::cout << "KVFExtractor::writeInfo(): pvHandle is not valid..." << std::endl;
              }

              vector< RefCountedKinematicParticle > jpsi1_children = vertexFitTree->finalStateParticles();
              if (jpsi1_children.size() != 2) {
                std::cout << " Warning Jpsi1 children size not equal to 2..." << std::endl;
              } else {
                TLorentzVector p_mu;
                // Order is : x,y,z,px,py,pz,m
                AlgebraicVector7 par1 = jpsi1_children[0]->currentState().kinematicParameters().vector();
                double e1 = jpsi1_children[0]->currentState().kinematicParameters().energy();
                new((*m_jpsikvf_mu1_lorentzvector)[nJpsi-1]) TLorentzVector(par1(3),par1(4),par1(5),e1);
                p_mu.SetPxPyPzE(par1(3),par1(4),par1(5),e1);
                if (m_isMC) {
                  for (auto& it: m_jpsikvf_mu1_muon_scaleFactors) {
                    std::pair<ScaleFactorService::WorkingPoint, ScaleFactorService::WorkingPoint> workingPoints = ScaleFactorService::getWorkingPointFromName(it.first);
                    it.second.push_back(m_scaleFactorService->getMuonScaleFactor(workingPoints.first, workingPoints.second, p_mu.Pt(), p_mu.Eta()));
                  }
                }

                AlgebraicVector7 par2 = jpsi1_children[1]->currentState().kinematicParameters().vector();
                double e2 = jpsi1_children[1]->currentState().kinematicParameters().energy();
                new((*m_jpsikvf_mu2_lorentzvector)[nJpsi-1]) TLorentzVector(par2(3),par2(4),par2(5),e2);
                p_mu.SetPxPyPzE(par2(3),par2(4),par2(5),e2);
                if (m_isMC) {
                  for (auto& it: m_jpsikvf_mu2_muon_scaleFactors) {
                    std::pair<ScaleFactorService::WorkingPoint, ScaleFactorService::WorkingPoint> workingPoints = ScaleFactorService::getWorkingPointFromName(it.first);
                    it.second.push_back(m_scaleFactorService->getMuonScaleFactor(workingPoints.first, workingPoints.second, p_mu.Pt(), p_mu.Eta()));
                  }
                }
              }

            } 
            else 
              std::cout << "J/psi decay vertex Not valid" << std::endl;

          } // vertex is not valid
        } // mass condition

      } // end 2nd PF loop
    } // end 1st PF loop

    m_jpsi_size = nJpsi;

    // end of J/psi stuff

    sort(myPFs.begin(), myPFs.end(), mSorterPFs);  
    sort(myKPis.begin(), myKPis.end(), mSorterPFs);  

    if (hasNonIsoMu) {
      ++nMuJet;

      if (unfold) {
        for (unsigned int j = 0; j < (unsigned int)myPFs2Unfold.size(); j++) {
          ++nUnfoldTr;
          m_mujet_unfold_indmujet[nUnfoldTr-1] = nMuJet-1;
          TLorentzVector recoP, genP;
          recoP.SetPxPyPzE(myPFs2Unfold[j].px(), myPFs2Unfold[j].py(), myPFs2Unfold[j].pz(), myPFs2Unfold[j].energy());
          double dRmin = 200.;
          for (unsigned int k = 0; k < (unsigned int)genParticles->size(); k++) {
            if (((*genParticles)[k]).px() < 1e-6 && ((*genParticles)[k]).py() < 1e-6) continue;
            TLorentzVector genP_int;
            genP_int.SetPxPyPzE(((*genParticles)[k]).px(), ((*genParticles)[k]).py(), ((*genParticles)[k]).pz(), ((*genParticles)[k]).energy());
            if (genP_int.DeltaR(recoP) < dRmin) {
              dRmin = genP_int.DeltaR(recoP);
              genP.SetPxPyPzE(((*genParticles)[k]).px(), ((*genParticles)[k]).py(), ((*genParticles)[k]).pz(), ((*genParticles)[k]).energy());
            }
          }
          if (fabs(myPFs2Unfold[j].pdgId()) != 13 ) {
            m_mujet_unfold_tr_recopt[nUnfoldTr-1] = recoP.Pt();
            m_mujet_unfold_tr_recoeta[nUnfoldTr-1] = recoP.Eta();
            m_mujet_unfold_tr_genpt[nUnfoldTr-1] = genP.Pt();
            m_mujet_unfold_tr_geneta[nUnfoldTr-1] = genP.Eta();
            m_mujet_unfold_tr_dr[nUnfoldTr-1] = dRmin;
          }
          else {
            m_mujet_unfold_mu_recopt[nUnfoldTr-1] = recoP.Pt();
            m_mujet_unfold_mu_recoeta[nUnfoldTr-1] = recoP.Eta();
            m_mujet_unfold_mu_genpt[nUnfoldTr-1] = genP.Pt();
            m_mujet_unfold_mu_geneta[nUnfoldTr-1] = genP.Eta();
            m_mujet_unfold_mu_dr[nUnfoldTr-1] = dRmin;
          }
        }
      }

      m_mujet_jet_btag_CSV[nMuJet-1] = p_jets.at(i).bDiscriminator("combinedSecondaryVertexBJetTags");  
      new((*m_mujet_jet_lorentzvector)[nMuJet-1]) TLorentzVector((p_jets.at(i)).px(),(p_jets.at(i)).py(),(p_jets.at(i)).pz(),(p_jets.at(i)).energy()); 
      if (m_isMC) {
        int mcFlavor = abs(p_jets.at(i).partonFlavour());
        ScaleFactorService::Flavor flavor = ScaleFactorService::B;
        if (mcFlavor == 4) {
          flavor = ScaleFactorService::C;
        } else if ((mcFlavor <= 3) || (mcFlavor == 21)) {
          // If mcFlavor == 0, assume it's a light jet
          flavor = ScaleFactorService::LIGHT;
        }
        m_mujet_jet_scaleFactors.push_back(m_scaleFactorService->getBTaggingScaleFactor(flavor, p_jets.at(i).et(), p_jets.at(i).eta()));
      }

      m_mujet_ntr[nMuJet-1] = myPFs.size();
      m_mujet_sump[nMuJet-1] = SumP;
      m_mujet_sumpt[nMuJet-1] = SumPt;
      m_mujet_sumvecp[nMuJet-1] = (float)SumVecP.P(); 
      bool nomuplus = true;
      bool nomuminus = true;
      for (unsigned int j = 0; j < (unsigned int)myPFs.size(); j++) {
        if (myPFs[j].pdgId() == 13 && nomuplus) {
          nomuplus = false;
          new((*m_mujet_nonisomuplus_lorentzvector)[nMuJet-1]) TLorentzVector(myPFs[j].px(),myPFs[j].py(),myPFs[j].pz(),myPFs[j].energy());
          m_mujet_nonisomuplus_pdgid[nMuJet-1] = myPFs[j].pdgId();
          if (m_isMC) {
            for (auto& it: m_mujet_nonisomuplus_scaleFactors) {
              std::pair<ScaleFactorService::WorkingPoint, ScaleFactorService::WorkingPoint> workingPoints = ScaleFactorService::getWorkingPointFromName(it.first);
              it.second.push_back(m_scaleFactorService->getMuonScaleFactor(workingPoints.first, workingPoints.second, myPFs[j].pt(), myPFs[j].eta()));
            }
          }
        }
        if (myPFs[j].pdgId() == -13 && nomuminus) {
          nomuminus = false;
          new((*m_mujet_nonisomuminus_lorentzvector)[nMuJet-1]) TLorentzVector(myPFs[j].px(),myPFs[j].py(),myPFs[j].pz(),myPFs[j].energy());
          m_mujet_nonisomuminus_pdgid[nMuJet-1] = myPFs[j].pdgId();
          if (m_isMC) {
            for (auto& it: m_mujet_nonisomuminus_scaleFactors) {
              std::pair<ScaleFactorService::WorkingPoint, ScaleFactorService::WorkingPoint> workingPoints = ScaleFactorService::getWorkingPointFromName(it.first);
              it.second.push_back(m_scaleFactorService->getMuonScaleFactor(workingPoints.first, workingPoints.second, myPFs[j].pt(), myPFs[j].eta()));
            }
          }
        }
        if (m_mujet_nonisomuplus_pdgid[nMuJet-1] != 0 && m_mujet_nonisomuminus_pdgid[nMuJet-1] != 0) break;
      }
      if (nomuplus) {
        new((*m_mujet_nonisomuplus_lorentzvector)[nMuJet-1]) TLorentzVector(0.,0.,0.,0.);
        if (m_isMC) {
          for (auto& it: m_mujet_nonisomuplus_scaleFactors) {
            it.second.push_back(ScaleFactor(0,0,0));
          }
        }
      }
      if (nomuminus) {
        new((*m_mujet_nonisomuminus_lorentzvector)[nMuJet-1]) TLorentzVector(0.,0.,0.,0.);    
        if (m_isMC) {
          for (auto& it: m_mujet_nonisomuminus_scaleFactors) {
            it.second.push_back(ScaleFactor(0,0,0));
          }
        }      
      }
      for (unsigned int j = 0; j < std::min((unsigned int)myPFs.size(), (unsigned int)m_tr_MAX); j++) {
        ++nTr;
        new((*m_mujet_tr_lorentzvector)[nTr-1]) TLorentzVector(myPFs[j].px(),myPFs[j].py(),myPFs[j].pz(),myPFs[j].energy());
        m_mujet_tr_pdgid[nTr-1] = myPFs[j].pdgId();
        m_mujet_tr_indmujet[nTr-1] = nMuJet-1;
      }
      for (unsigned int j = std::min((unsigned int)myPFs.size(), (unsigned int)m_tr_MAX); j < (unsigned int)m_tr_MAX; j++) {
        ++nTr;
        new((*m_mujet_tr_lorentzvector)[nTr-1]) TLorentzVector(0., 0., 0., 0.);
        m_mujet_tr_pdgid[nTr-1] = 0;
        m_mujet_tr_indmujet[nTr-1] = nMuJet-1;
      }


      // Reconstruct the D0
      //----------------------
      ParticleMass kaon_mass  = 0.493677;
      float        kaon_sigma = 0.000001;

      ParticleMass pion_mass  = 0.13957018;
      float        pion_sigma = 0.00000001;

      // ParticleMass d0_mass  = 1.86484;
      // float        d0_sigma = 0.00014;

      // To transform Track to TransientTrack, first need to get the builder:
      edm::ESHandle<TransientTrackBuilder> theKPiB;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theKPiB);  

      int nd0 = 0;

      for(unsigned int j = 0; j < std::min((unsigned int)myKPis.size(), (unsigned int)m_nTrD0Max); ++j) {

        for(unsigned int k = 0; k < std::min((unsigned int)myKPis.size(), (unsigned int)m_nTrD0Max); ++k) {

          // PF particles should be different and of OS
          if (j == k) continue;
          if (myKPis[j].charge() + myKPis[k].charge() != 0) continue;

          ++nd0;
          ++nD0;
          m_mujet_d0kvf_indmujet[nD0-1] = nMuJet-1;

          double eD0  = myKPis[j].energy()+myKPis[k].energy();
          double pxD0 = myKPis[j].px()+myKPis[k].px();
          double pyD0 = myKPis[j].py()+myKPis[k].py();
          double pzD0 = myKPis[j].pz()+myKPis[k].pz();
          double mD0 = pow(eD0,2)-pow(pxD0,2)-pow(pyD0,2)-pow(pzD0,2);
          if (mD0 > 0.) mD0 = sqrt(mD0);
          else mD0 = 0.;

          // save the D0 from 2 PF particles (which can be different from D0 from 2 tracks)
          new((*m_mujet_d0pf_lorentzvector)[nD0-1]) TLorentzVector(pxD0, pyD0, pzD0, eD0);

          // Make the Transient tracks
          reco::TransientTrack tr1 = (*theKPiB).build(myKPis[j].trackRef());// j one is kaon
          reco::TransientTrack tr2 = (*theKPiB).build(myKPis[k].trackRef());// k one is pion

          // A complicated Kalman fit :
          //---------------------------

          //Creating a KinematicParticleFactory
          KinematicParticleFactoryFromTransientTrack pFactory;

          //initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered
          float chi = 0.;
          float ndf = 0.;

          //making particles
          std::vector<RefCountedKinematicParticle> kpiParticles;
          kpiParticles.push_back(pFactory.particle (tr1,kaon_mass,chi,ndf,kaon_sigma));
          kpiParticles.push_back(pFactory.particle (tr2,pion_mass,chi,ndf,pion_sigma));

          /* Example of a simple vertex fit, without other constraints
           * The reconstructed decay tree is a result of the kinematic fit
           * The KinematicParticleVertexFitter fits the final state particles to their vertex and
           * reconstructs the decayed state
           */

          // creating the vertex fitter
          KinematicParticleVertexFitter fitter;
          // reconstructing a J/Psi decay
          RefCountedKinematicTree vertexFitTree = fitter.fit(kpiParticles);

          if (!vertexFitTree->isValid()) {
            //std::cout <<"D0 vertexTree is invalid. Fit failed." << std::endl;

            // Need to fill empty quantities : OK for tables, but need to create empty TLorentzVector
            new((*m_mujet_d0kvf_lorentzvector)[nD0-1]) TLorentzVector(0.,0.,0.,0.);
            new((*m_mujet_d0kvf_pion_lorentzvector)[nD0-1]) TLorentzVector(0.,0.,0.,0.);
            new((*m_mujet_d0kvf_kaon_lorentzvector)[nD0-1]) TLorentzVector(0.,0.,0.,0.);

            continue;

          } else {
            //accessing the tree components, move pointer to top
            vertexFitTree->movePointerToTheTop();

            //We are now at the top of the decay tree getting the jpsi reconstructed KinematicPartlcle
            RefCountedKinematicParticle d0 = vertexFitTree->currentParticle();
            AlgebraicVector7 par0 = d0->currentState().kinematicParameters().vector();
            double e0 = d0->currentState().kinematicParameters().energy();
            new((*m_mujet_d0kvf_lorentzvector)[nD0-1]) TLorentzVector(par0(3),par0(4),par0(5),e0);

            RefCountedKinematicVertex d0_vertex = vertexFitTree->currentDecayVertex();
            if ( d0_vertex->vertexIsValid()) {

              m_mujet_d0kvf_vx[nD0-1] = d0_vertex->position().x();
              m_mujet_d0kvf_vy[nD0-1] = d0_vertex->position().y();
              m_mujet_d0kvf_vz[nD0-1] = d0_vertex->position().z();
              m_mujet_d0kvf_vtxvalid[nD0-1] = true;
              m_mujet_d0kvf_vtxchi2[nD0-1] = d0_vertex->chiSquared();
              m_mujet_d0kvf_ndf[nD0-1] = d0_vertex->degreesOfFreedom();

              // Compute the distance between the PV and the D0 vertex :
              //----------------------------------------------------------

              const reco::VertexCollection vtx = *(pvHandle.product());

              GlobalPoint svPos    = d0_vertex->position();
              GlobalError svPosErr = d0_vertex->error();

              // If the  PV does not exist, compute distance wrt to the detector center (0,0,0)

              double sigmax = 0.;
              double sigmay = 0.;
              double sigmaz = 0.;
              if (vtx.size() > 0) {
                sigmax = sqrt(vtx[0].xError()*vtx[0].xError() + svPosErr.cxx()*svPosErr.cxx());
                sigmay = sqrt(vtx[0].yError()*vtx[0].yError() + svPosErr.cyy()*svPosErr.cyy());
                sigmaz = sqrt(vtx[0].zError()*vtx[0].zError() + svPosErr.czz()*svPosErr.czz());
              }  else {
                sigmax = sqrt(svPosErr.cxx()*svPosErr.cxx());
                sigmay = sqrt(svPosErr.cyy()*svPosErr.cyy());
                sigmaz = sqrt(svPosErr.czz()*svPosErr.czz());
              }

              double px  = par0(3);
              double py  = par0(4);
              double pz  = par0(5);
              double nrj = e0;
              double m = sqrt(nrj*nrj - px*px - py*py - pz*pz);

              double interx = pow((px/m)/sigmax, 2.);
              double intery = pow((py/m)/sigmay, 2.);
              double interz = pow((pz/m)/sigmaz, 2.);

              m_mujet_d0kvf_sigmaL3D[nD0-1] = pow( interx + intery + interz , -0.5);
              //std::cout << "sigmaL3D = " << m_mujet_d0kvf_sigmaL3D[nD0-1] << std::endl;

              double part1 = 0.;
              double part2 = 0.;
              double part3 = 0.;
              if (vtx.size() >0) {
                part1 = (px/m)*pow(m_mujet_d0kvf_sigmaL3D[nD0-1]/sigmax,2.)*(svPos.x() - vtx[0].x());
                part2 = (py/m)*pow(m_mujet_d0kvf_sigmaL3D[nD0-1]/sigmay,2.)*(svPos.y() - vtx[0].y());
                part3 = (pz/m)*pow(m_mujet_d0kvf_sigmaL3D[nD0-1]/sigmaz,2.)*(svPos.z() - vtx[0].z());
              }  else {
                part1 = (px/m)*pow(m_mujet_d0kvf_sigmaL3D[nD0-1]/sigmax,2.)*svPos.x();
                part2 = (py/m)*pow(m_mujet_d0kvf_sigmaL3D[nD0-1]/sigmay,2.)*svPos.y();
                part3 = (pz/m)*pow(m_mujet_d0kvf_sigmaL3D[nD0-1]/sigmaz,2.)*svPos.z();
              }

              m_mujet_d0kvf_L3D[nD0-1] = fabs(part1 + part2 + part3);
              //std::cout << "L3D = " << m_mujet_d0kvf_L3D[nD0-1] << std::endl;  

              m_mujet_d0kvf_L3DoverSigmaL3D[nD0-1] = m_mujet_d0kvf_L3D[nD0-1]/m_mujet_d0kvf_sigmaL3D[nD0-1];
              //std::cout << "(L/sigma)3D = " << m_mujet_d0kvf_L3DoverSigmaL3D[nD0-1] << std::endl;    

              if (!pvHandle.isValid()) {
                std::cout << "KVFExtractor::writeInfo(): pvHandle is not valid..." << std::endl;
              }

              vector<RefCountedKinematicParticle> d0_children = vertexFitTree->finalStateParticles();
              if (d0_children.size() != 2) {
                std::cout << " Warning D0 children size not equal to 2..." << std::endl;
              } else {
                // Order is : x,y,z,px,py,pz,m
                AlgebraicVector7 par1 = d0_children[0]->currentState().kinematicParameters().vector();
                double e1 = d0_children[0]->currentState().kinematicParameters().energy();
                new((*m_mujet_d0kvf_kaon_lorentzvector)[nD0-1]) TLorentzVector(par1(3),par1(4),par1(5),e1);
                m_mujet_d0kvf_kaon_pdgid[nD0-1] = myKPis[j].charge()*321; 

                AlgebraicVector7 par2 = d0_children[1]->currentState().kinematicParameters().vector();
                double e2 = d0_children[1]->currentState().kinematicParameters().energy();
                new((*m_mujet_d0kvf_pion_lorentzvector)[nD0-1]) TLorentzVector(par2(3),par2(4),par2(5),e2);
                m_mujet_d0kvf_pion_pdgid[nD0-1] = myKPis[k].charge()*211; 
              }

            } 
            else 
              std::cout << "D0 decay vertex Not valid" << std::endl;

          } // vertex is not valid

        } // end 2nd PF loop
      } // end 1st PF loop

      m_mujet_nd0[nMuJet-1] = nd0;

      // end of D0 stuff

    } // mu tagged jet  

    m_size++;
  }  // jet loop

  m_mujet_size = nMuJet;
  m_mujet_tr_size = nTr;
  m_mujet_unfold_tr_size = nUnfoldTr;
  m_mujet_d0_size = nD0;

  // end of mu tagged jet stuff

  fillTree();
}

//
// Method getting the info from an input file
//

void KVFExtractor::getInfo(int ievt) 
{
  if (m_tree_jpsi)
    m_tree_jpsi->GetEntry(ievt);

}


// Method initializing everything (to do for each event)

void KVFExtractor::reset()
{
  m_size = 0;

  // Jpsi tree
  m_jpsi_jet_lorentzvector->Clear();
  m_jpsi_jet_scaleFactors.clear();
  m_jpsipf_lorentzvector->Clear();
  m_jpsikvf_lorentzvector->Clear();
  m_jpsikvf_mu1_lorentzvector->Clear();
  m_jpsikvf_mu2_lorentzvector->Clear();

  for (auto& it: m_jpsikvf_mu1_muon_scaleFactors) {
    it.second.clear();
  }
  for (auto& it: m_jpsikvf_mu2_muon_scaleFactors) {
    it.second.clear();
  }

  m_jpsi_size = 0;

  for (int i = 0; i < m_jpsi_MAX; ++i) {
    m_jpsi_indjet[i] = 0;
    m_jpsi_jet_btag_CSV[i] = 0;
    m_jpsi_indpf1[i] = 0;
    m_jpsi_indpf2[i] = 0;

    m_jpsikvf_vx[i] = 0.;
    m_jpsikvf_vy[i] = 0.;
    m_jpsikvf_vz[i] = 0.;
    m_jpsikvf_vtxvalid[i] = false;
    m_jpsikvf_vtxchi2[i] = 0.;
    m_jpsikvf_ndf[i] = 0.;

    m_jpsikvf_L3D[i] = 0.;
    m_jpsikvf_sigmaL3D[i] = 0.;
    m_jpsikvf_L3DoverSigmaL3D[i] = 0.;
  }

  // mujet tree

  m_mujet_jet_lorentzvector->Clear();
  m_mujet_jet_scaleFactors.clear();
  m_mujet_nonisomuplus_lorentzvector->Clear();
  m_mujet_nonisomuminus_lorentzvector->Clear();
  m_mujet_tr_lorentzvector->Clear(); 
  m_mujet_d0pf_lorentzvector->Clear();
  m_mujet_d0kvf_lorentzvector->Clear();
  m_mujet_d0kvf_pion_lorentzvector->Clear();
  m_mujet_d0kvf_kaon_lorentzvector->Clear();

  for (auto& it: m_mujet_nonisomuplus_scaleFactors) {
    it.second.clear();
  }
  for (auto& it: m_mujet_nonisomuminus_scaleFactors) {
    it.second.clear();
  }

  m_mujet_size = 0;
  m_mujet_tr_size = 0;
  m_mujet_unfold_tr_size = 0;
  m_mujet_d0_size = 0;

  for (int i = 0; i < m_mujet_MAX; ++i) {
    m_mujet_jet_btag_CSV[i] = 0;
    m_mujet_nonisomuplus_pdgid[i] = 0;
    m_mujet_nonisomuminus_pdgid[i] = 0;
    m_mujet_ntr[i] = 0;
    m_mujet_sump[i] = 0;
    m_mujet_sumpt[i] = 0;
    m_mujet_sumvecp[i] = 0;
    m_mujet_nd0[i] = 0;
  }
  for (int i = 0; i < m_Tr_MAX; i++) {
    m_mujet_tr_indmujet[i] = 0;
    m_mujet_tr_pdgid[i] = 0;
  } 
  for (int i = 0; i < m_D0_MAX; ++i) {
    m_mujet_d0kvf_indmujet[i] = 0;
    m_mujet_d0kvf_pion_pdgid[i] = 0; 
    m_mujet_d0kvf_kaon_pdgid[i] = 0; 
    m_mujet_d0kvf_vx[i] = 0;  
    m_mujet_d0kvf_vy[i] = 0;  
    m_mujet_d0kvf_vz[i] = 0;
    m_mujet_d0kvf_vtxvalid[i] = 0;
    m_mujet_d0kvf_vtxchi2[i] = 0;
    m_mujet_d0kvf_ndf[i] = 0;  
    m_mujet_d0kvf_L3D[i] = 0;  
    m_mujet_d0kvf_sigmaL3D[i] = 0;  
    m_mujet_d0kvf_L3DoverSigmaL3D[i] = 0;  
  }
  for (int i = 0; i < m_unfold_Tr_MAX; ++i) {
    m_mujet_unfold_indmujet[i] = 0;
    m_mujet_unfold_tr_recopt[i] = 0;
    m_mujet_unfold_tr_recoeta[i] = 0;
    m_mujet_unfold_tr_genpt[i] = 0;
    m_mujet_unfold_tr_geneta[i] = 0;
    m_mujet_unfold_tr_dr[i] = 0;
    m_mujet_unfold_mu_recopt[i] = 0;
    m_mujet_unfold_mu_recoeta[i] = 0;
    m_mujet_unfold_mu_genpt[i] = 0;
    m_mujet_unfold_mu_geneta[i] = 0;
    m_mujet_unfold_mu_dr[i] = 0;
  }
}


void KVFExtractor::fillTree()
{
  if (m_tree_jpsi)
    m_tree_jpsi->Fill(); 
  if (m_tree_mujet)
    m_tree_mujet->Fill(); 

}

void KVFExtractor::correctJets(pat::JetCollection& jets, const edm::Event& iEvent, const edm::EventSetup& iSetup) {

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Recompute jet energy corrections..." << std::endl;
#endif

  // Get Jet corrector
  const JetCorrector* globalTagCorrector = nullptr;

  if (mUseGlobalTagForJEC) {
    globalTagCorrector = JetCorrector::getJetCorrector(mJetCorrectorLabel, iSetup);
  }

  edm::Handle<reco::VertexCollection>  vertexHandle;
  iEvent.getByLabel("goodOfflinePrimaryVertices", vertexHandle);
  reco::VertexCollection vertices = *vertexHandle;

  edm::Handle<double> rhos;
  iEvent.getByLabel(edm::InputTag("kt6PFJets", "rho", "RECO"), rhos);
  double rho = *rhos;

  // Correct jets
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    pat::Jet& jet = *it;

    pat::Jet rawJet = jet.correctedJet("Uncorrected");
    jet.addUserData("rawJet", rawJet, true); // Store raw jet inside the jet
    const pat::Jet L1Jet  = jet.correctedJet("L1FastJet");
    jet.addUserData("L1Jet", L1Jet, true); // Embed L1 corrected jet for TypeI correction

#if DEBUG
    std::cout << "---" << std::endl;
    std::cout << "Pt: " << jet.pt() << std::endl;
#endif

    double toRaw = jet.jecFactor("Uncorrected");
    jet.setP4(jet.p4() * toRaw); // jet is now a raw jet
#if DEBUG
    std::cout << "True raw pt: " << rawJet.pt() << std::endl;
    std::cout << "Raw pt: " << jet.pt() << std::endl;
#endif

    double corrections = 0.;
    if (mUseGlobalTagForJEC) {
      corrections = globalTagCorrector->correction(jet, iEvent, iSetup);
    } else {
      mTxtCorrector->setJetEta(jet.eta());
      mTxtCorrector->setJetPt(jet.pt());
      mTxtCorrector->setRho(rho);
      mTxtCorrector->setJetA(jet.jetArea());
      mTxtCorrector->setNPV(vertices.size());
      corrections = mTxtCorrector->getCorrection();
    }

    jet.scaleEnergy(corrections);
#if DEBUG
    std::cout << "Corrected pt: " << jet.pt() << std::endl;
#endif
  }

  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorterJets);
}

//from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution

double KVFExtractor::getResCorrFactor(const pat::Jet& jet) {

  double factor = 0.;
  double error  = 0.;

  if (fabs(jet.eta()) > 0. && fabs(jet.eta()) <= 0.5) {
    factor = 1.052;
    error = (mJERSign == 1) ? 0.062 : 0.061;
  } else if (fabs(jet.eta()) > 0.5 && fabs(jet.eta()) <= 1.1) {
    factor = 1.057;
    error = (mJERSign == 1) ? 0.056 : 0.055;
  } else if (fabs(jet.eta()) > 1.1 && fabs(jet.eta()) <= 1.7) {
    factor = 1.096;
    error = (mJERSign == 1) ? 0.063 : 0.062;
  } else if (fabs(jet.eta()) > 1.7 && fabs(jet.eta()) <= 2.3) {
    factor = 1.134;
    error = (mJERSign == 1) ? 0.087 : 0.085;
  } else if (fabs(jet.eta()) > 2.3 && fabs(jet.eta()) <= 5.0) {
    factor = 1.288;
    error = (mJERSign == 1) ? 0.155 : 0.153;
  }

  return factor + mJERSign * error;
}

void KVFExtractor::correctJetsResolution(pat::JetCollection& jets) {

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Doing jet resolution smearing" << std::endl;
#endif

  // Correct jets
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    pat::Jet& jet = *it;

    if (jet.pt() > 10) {

#if DEBUG
      std::cout << "---" << std::endl;
      std::cout << "Pt: " << jet.pt() << std::endl;
#endif

      // resolution corection as in https://twiki.cern.ch/twiki/bin/view/CMS/TWikiTopRefSyst#Jet_energy_resolution

      if (jet.genJet() == NULL)
        continue;

      double genjet_pt = jet.genJet()->pt();
      double jet_pt = jet.pt();
      double rescorr = getResCorrFactor(jet) - 1;
      double deltapt = (jet_pt - genjet_pt) * rescorr; 
      double scalefac = (jet_pt + deltapt) / jet_pt;
      if (scalefac <= 0)
        continue;

      jet.scaleEnergy(scalefac);

#if DEBUG
      std::cout << "Corrected pt: " << jet.pt() << std::endl;
#endif

    }
  }

  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorterJets);
}


//from JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi, check updates
//updated to MET systematic shift corrections to 2012 ABCD ReReco data + new Summer'13 JEC
double KVFExtractor::getSysShifCorrFactorX(const int Nvtx){
  if (m_isMC) return -(+1.62861e-01 - 2.38517e-02*Nvtx);
  else        return -(+4.83642e-02 + 2.48870e-01*Nvtx);
}

double KVFExtractor::getSysShifCorrFactorY(const int Nvtx){
  if (m_isMC) return -(+3.60860e-01 - 1.30335e-01*Nvtx);
  else        return -(-1.50135e-01 - 8.27917e-02*Nvtx);
}


void KVFExtractor::doJESSystematics(pat::JetCollection& jets) {
#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "JES systematic" << std::endl;
#endif

  // Correct jets
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    pat::Jet& jet = *it;

    jecUncertainty->setJetEta(jet.eta());
    jecUncertainty->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt

    double uncertainty = (mJESSign == 1) ? fabs(jecUncertainty->getUncertainty(true)) : fabs(jecUncertainty->getUncertainty(false));
    double signedCorrection = mJESSign * uncertainty;

    double scaleFactor = 1. + signedCorrection;

#if DEBUG
    std::cout << "---" << std::endl;
    std::cout << "Pt before JES uncertainty: " << jet.pt() << std::endl;
#endif

    jet.scaleEnergy(scaleFactor);

#if DEBUG
    std::cout << "Pt after JES uncertainty: " << jet.pt() << std::endl;
#endif

  }

  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorterJets);
}

void KVFExtractor::extractRawJets(pat::JetCollection& jets) {

  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it) {
    pat::Jet& jet = *it;

    const pat::Jet rawJet = jet.correctedJet("Uncorrected");
    jet.addUserData("rawJet", rawJet, true);
    const pat::Jet L1Jet  = jet.correctedJet("L1FastJet");
    jet.addUserData("L1Jet", L1Jet, true); // Embed L1 corrected jet for TypeI correction
  }

}

