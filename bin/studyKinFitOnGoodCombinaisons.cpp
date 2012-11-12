/**
 * Study KinFit on mtt mass, only on good combinaisons (associable AND matched)
 */


//#include <Riostream.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <algorithm>
#include <string>

#include <ios>

#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <TH1F.h>
#include <TF1.h>
#include <TROOT.h>
#include <TCut.h>
#include <TProfile.h>
#include <TError.h>

#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/PythonParameterSet/interface/MakeParameterSets.h>
#include <TopQuarkAnalysis/TopKinFitter/interface/CovarianceMatrix.h>

#include "FWCore/MessageService/interface/MessageServicePresence.h"

#include "KinFitter.h"
#include "GaussianProfile.h"

#include "tclap/CmdLine.h"

// Variables
static const int 	m_jets_MAX       = 200;

TClonesArray* m_jet_lorentzvector = new TClonesArray("TLorentzVector");
int  n_jets;
int  m_jet_MCIndex[m_jets_MAX];

float MC_mtt;

int genLeptonicBIndex, genHadronicBIndex, genFirstJetIndex, genSecondJetIndex, genLeptonIndex, genNeutrinoIndex;

const int nBins = 15;
const double bins[] = {340, 360, 380, 400, 420, 460, 500, 550, 600, 650, 750, 850, 950, 1050, 1200, 1400, 1600};

bool CSV_MODE = true;
uint64_t ENTRIES = 0;

double binomialError(double efficiency, int n) {
  return sqrt((efficiency * (1 - efficiency)) / (n));
}

void loadChain(const std::vector<std::string>& inputFiles, TChain*& mc, TChain*& event, TChain*& jets, TChain*& MET, TChain*& muons, TChain*& electrons, TChain*& mtt) {
  std::cout << "Opening files..." << std::endl;

  mc = new TChain("MC");
  event = new TChain("event");
  jets = new TChain("jet_PF");
  MET = new TChain("MET_PF");
  muons = new TChain("muon_PF");
  electrons = new TChain("electron_PF");
  mtt = new TChain("Mtt");

  for (const std::string& file: inputFiles) {
    mc->Add(file.c_str());
    event->Add(file.c_str());
    jets->Add(file.c_str());
    MET->Add(file.c_str());
    muons->Add(file.c_str());
    electrons->Add(file.c_str());
    mtt->Add(file.c_str());
  }

  std::cout << "... done." << std::endl;
}

void SetBranchAddress(TChain* chain, const char* branchName, void* address) {
  chain->SetBranchStatus(branchName, 1);
  chain->SetBranchAddress(branchName, address, NULL);
}

bool isMatched(int index, int jetsMcIndex[], TClonesArray* p4s, int nJets) {
  for (int i = 0; i < nJets; i++) {
    if (jetsMcIndex[i] == index) {
      TLorentzVector* p4 = static_cast<TLorentzVector*>((*p4s)[i]);
      return p4->Pt() > 30 && fabs(p4->Eta()) < 2.4;
    }
  }

  return false;
}

void process(const std::vector<std::string>& inputFiles, const std::string& outputFile, const edm::ParameterSet& kinFitPS) {

  //TH1::SetDefaultSumw2(true);
  gROOT->SetBatch(true);

  KinFitter kinFitter(kinFitPS);
  TtSemiLepKinFitter& fitter = kinFitter.getFitter();
  fitter.setVerbosity(0);

  TTree* tree = new TTree("tree", "tree");

  float mtt_reco;
  tree->Branch("mtt", &mtt_reco, "mtt/F");

  float mtt_reco_after_kf;
  tree->Branch("mtt_after_kf", &mtt_reco_after_kf, "mtt_after_kf/F");

  float firstJet_pt, firstJet_eta, firstJet_phi, firstJet_et, firstJet_px, firstJet_py, firstJet_pz;
  tree->Branch("firstJet_pt", &firstJet_pt, "firstJet_pt/F");
  tree->Branch("firstJet_eta", &firstJet_eta, "firstJet_eta/F");
  tree->Branch("firstJet_phi", &firstJet_phi, "firstJet_phi/F");
  tree->Branch("firstJet_et", &firstJet_et, "firstJet_et/F");
  tree->Branch("firstJet_px", &firstJet_px, "firstJet_px/F");
  tree->Branch("firstJet_py", &firstJet_py, "firstJet_py/F");
  tree->Branch("firstJet_pz", &firstJet_pz, "firstJet_pz/F");
  float firstJet_E;
  tree->Branch("firstJet_E", &firstJet_E, "firstJet_E/F");

  float secondJet_pt, secondJet_eta, secondJet_phi, secondJet_et, secondJet_px, secondJet_py, secondJet_pz;
  tree->Branch("secondJet_pt", &secondJet_pt, "secondJet_pt/F");
  tree->Branch("secondJet_eta", &secondJet_eta, "secondJet_eta/F");
  tree->Branch("secondJet_phi", &secondJet_phi, "secondJet_phi/F");
  tree->Branch("secondJet_et", &secondJet_et, "secondJet_et/F");
  tree->Branch("secondJet_px", &secondJet_px, "secondJet_px/F");
  tree->Branch("secondJet_py", &secondJet_py, "secondJet_py/F");
  tree->Branch("secondJet_pz", &secondJet_pz, "secondJet_pz/F");
  float secondJet_E;
  tree->Branch("secondJet_E", &secondJet_E, "secondJet_E/F");

  float hadronicB_pt, hadronicB_eta, hadronicB_phi, hadronicB_et, hadronicB_px, hadronicB_py, hadronicB_pz;
  tree->Branch("hadronicB_pt", &hadronicB_pt, "hadronicB_pt/F");
  tree->Branch("hadronicB_eta", &hadronicB_eta, "hadronicB_eta/F");
  tree->Branch("hadronicB_phi", &hadronicB_phi, "hadronicB_phi/F");
  tree->Branch("hadronicB_et", &hadronicB_et, "hadronicB_et/F");
  tree->Branch("hadronicB_px", &hadronicB_px, "hadronicB_px/F");
  tree->Branch("hadronicB_py", &hadronicB_py, "hadronicB_py/F");
  tree->Branch("hadronicB_pz", &hadronicB_pz, "hadronicB_pz/F");
  float hadronicB_E;
  tree->Branch("hadronicB_E", &hadronicB_E, "hadronicB_E/F");

  float leptonicB_pt, leptonicB_eta, leptonicB_phi, leptonicB_et, leptonicB_px, leptonicB_py, leptonicB_pz;
  tree->Branch("leptonicB_pt", &leptonicB_pt, "leptonicB_pt/F");
  tree->Branch("leptonicB_eta", &leptonicB_eta, "leptonicB_eta/F");
  tree->Branch("leptonicB_phi", &leptonicB_phi, "leptonicB_phi/F");
  tree->Branch("leptonicB_et", &leptonicB_et, "leptonicB_et/F");
  tree->Branch("leptonicB_px", &leptonicB_px, "leptonicB_px/F");
  tree->Branch("leptonicB_py", &leptonicB_py, "leptonicB_py/F");
  tree->Branch("leptonicB_pz", &leptonicB_pz, "leptonicB_pz/F");
  float leptonicB_E;
  tree->Branch("leptonicB_E", &leptonicB_E, "leptonicB_E/F");

  float neutrino_pt, neutrino_eta, neutrino_phi, neutrino_et, neutrino_px, neutrino_py, neutrino_pz;
  tree->Branch("neutrino_pt", &neutrino_pt, "neutrino_pt/F");
  tree->Branch("neutrino_eta", &neutrino_eta, "neutrino_eta/F");
  tree->Branch("neutrino_phi", &neutrino_phi, "neutrino_phi/F");
  tree->Branch("neutrino_et", &neutrino_et, "neutrino_et/F");
  tree->Branch("neutrino_px", &neutrino_px, "neutrino_px/F");
  tree->Branch("neutrino_py", &neutrino_py, "neutrino_py/F");
  tree->Branch("neutrino_pz", &neutrino_pz, "neutrino_pz/F");
  float neutrino_E;
  tree->Branch("neutrino_E", &neutrino_E, "neutrino_E/F");

  float lepton_pt, lepton_eta, lepton_phi, lepton_et, lepton_px, lepton_py, lepton_pz;
  tree->Branch("lepton_pt", &lepton_pt, "lepton_pt/F");
  tree->Branch("lepton_eta", &lepton_eta, "lepton_eta/F");
  tree->Branch("lepton_phi", &lepton_phi, "lepton_phi/F");
  tree->Branch("lepton_et", &lepton_et, "lepton_et/F");
  tree->Branch("lepton_px", &lepton_px, "lepton_px/F");
  tree->Branch("lepton_py", &lepton_py, "lepton_py/F");
  tree->Branch("lepton_pz", &lepton_pz, "lepton_pz/F");
  float lepton_E;
  tree->Branch("lepton_E", &lepton_E, "lepton_E/F");

  float hadronic_top_mass, leptonic_top_mass, hadronic_w_mass, leptonic_w_mass;
  tree->Branch("hadronic_top_mass", &hadronic_top_mass, "hadronic_top_mass/F");
  tree->Branch("leptonic_top_mass", &leptonic_top_mass, "leptonic_top_mass/F");
  tree->Branch("hadronic_w_mass", &hadronic_w_mass, "hadronic_w_mass/F");
  tree->Branch("leptonic_w_mass", &leptonic_w_mass, "leptonic_w_mass/F");

  float firstJet_after_kf_pt, firstJet_after_kf_eta, firstJet_after_kf_phi, firstJet_after_kf_et, firstJet_after_kf_px, firstJet_after_kf_py, firstJet_after_kf_pz;
  tree->Branch("firstJet_after_kf_pt", &firstJet_after_kf_pt, "firstJet_after_kf_pt/F");
  tree->Branch("firstJet_after_kf_eta", &firstJet_after_kf_eta, "firstJet_after_kf_eta/F");
  tree->Branch("firstJet_after_kf_phi", &firstJet_after_kf_phi, "firstJet_after_kf_phi/F");
  tree->Branch("firstJet_after_kf_et", &firstJet_after_kf_et, "firstJet_after_kf_et/F");
  tree->Branch("firstJet_after_kf_px", &firstJet_after_kf_px, "firstJet_after_kf_px/F");
  tree->Branch("firstJet_after_kf_py", &firstJet_after_kf_py, "firstJet_after_kf_py/F");
  tree->Branch("firstJet_after_kf_pz", &firstJet_after_kf_pz, "firstJet_after_kf_pz/F");
  float firstJet_after_kf_E;
  tree->Branch("firstJet_after_kf_E", &firstJet_after_kf_E, "firstJet_after_kf_E/F");

  float secondJet_after_kf_pt, secondJet_after_kf_eta, secondJet_after_kf_phi, secondJet_after_kf_et, secondJet_after_kf_px, secondJet_after_kf_py, secondJet_after_kf_pz;
  tree->Branch("secondJet_after_kf_pt", &secondJet_after_kf_pt, "secondJet_after_kf_pt/F");
  tree->Branch("secondJet_after_kf_eta", &secondJet_after_kf_eta, "secondJet_after_kf_eta/F");
  tree->Branch("secondJet_after_kf_phi", &secondJet_after_kf_phi, "secondJet_after_kf_phi/F");
  tree->Branch("secondJet_after_kf_et", &secondJet_after_kf_et, "secondJet_after_kf_et/F");
  tree->Branch("secondJet_after_kf_px", &secondJet_after_kf_px, "secondJet_after_kf_px/F");
  tree->Branch("secondJet_after_kf_py", &secondJet_after_kf_py, "secondJet_after_kf_py/F");
  tree->Branch("secondJet_after_kf_pz", &secondJet_after_kf_pz, "secondJet_after_kf_pz/F");
  float secondJet_after_kf_E;
  tree->Branch("secondJet_after_kf_E", &secondJet_after_kf_E, "secondJet_after_kf_E/F");

  float hadronicB_after_kf_pt, hadronicB_after_kf_eta, hadronicB_after_kf_phi, hadronicB_after_kf_et, hadronicB_after_kf_px, hadronicB_after_kf_py, hadronicB_after_kf_pz;
  tree->Branch("hadronicB_after_kf_pt", &hadronicB_after_kf_pt, "hadronicB_after_kf_pt/F");
  tree->Branch("hadronicB_after_kf_eta", &hadronicB_after_kf_eta, "hadronicB_after_kf_eta/F");
  tree->Branch("hadronicB_after_kf_phi", &hadronicB_after_kf_phi, "hadronicB_after_kf_phi/F");
  tree->Branch("hadronicB_after_kf_et", &hadronicB_after_kf_et, "hadronicB_after_kf_et/F");
  tree->Branch("hadronicB_after_kf_px", &hadronicB_after_kf_px, "hadronicB_after_kf_px/F");
  tree->Branch("hadronicB_after_kf_py", &hadronicB_after_kf_py, "hadronicB_after_kf_py/F");
  tree->Branch("hadronicB_after_kf_pz", &hadronicB_after_kf_pz, "hadronicB_after_kf_pz/F");
  float hadronicB_after_kf_E;
  tree->Branch("hadronicB_after_kf_E", &hadronicB_after_kf_E, "hadronicB_after_kf_E/F");

  float leptonicB_after_kf_pt, leptonicB_after_kf_eta, leptonicB_after_kf_phi, leptonicB_after_kf_et, leptonicB_after_kf_px, leptonicB_after_kf_py, leptonicB_after_kf_pz;
  tree->Branch("leptonicB_after_kf_pt", &leptonicB_after_kf_pt, "leptonicB_after_kf_pt/F");
  tree->Branch("leptonicB_after_kf_eta", &leptonicB_after_kf_eta, "leptonicB_after_kf_eta/F");
  tree->Branch("leptonicB_after_kf_phi", &leptonicB_after_kf_phi, "leptonicB_after_kf_phi/F");
  tree->Branch("leptonicB_after_kf_et", &leptonicB_after_kf_et, "leptonicB_after_kf_et/F");
  tree->Branch("leptonicB_after_kf_px", &leptonicB_after_kf_px, "leptonicB_after_kf_px/F");
  tree->Branch("leptonicB_after_kf_py", &leptonicB_after_kf_py, "leptonicB_after_kf_py/F");
  tree->Branch("leptonicB_after_kf_pz", &leptonicB_after_kf_pz, "leptonicB_after_kf_pz/F");
  float leptonicB_after_kf_E;
  tree->Branch("leptonicB_after_kf_E", &leptonicB_after_kf_E, "leptonicB_after_kf_E/F");

  float neutrino_after_kf_pt, neutrino_after_kf_eta, neutrino_after_kf_phi, neutrino_after_kf_et, neutrino_after_kf_px, neutrino_after_kf_py, neutrino_after_kf_pz;
  tree->Branch("neutrino_after_kf_pt", &neutrino_after_kf_pt, "neutrino_after_kf_pt/F");
  tree->Branch("neutrino_after_kf_eta", &neutrino_after_kf_eta, "neutrino_after_kf_eta/F");
  tree->Branch("neutrino_after_kf_phi", &neutrino_after_kf_phi, "neutrino_after_kf_phi/F");
  tree->Branch("neutrino_after_kf_et", &neutrino_after_kf_et, "neutrino_after_kf_et/F");
  tree->Branch("neutrino_after_kf_px", &neutrino_after_kf_px, "neutrino_after_kf_px/F");
  tree->Branch("neutrino_after_kf_py", &neutrino_after_kf_py, "neutrino_after_kf_py/F");
  tree->Branch("neutrino_after_kf_pz", &neutrino_after_kf_pz, "neutrino_after_kf_pz/F");
  float neutrino_after_kf_E;
  tree->Branch("neutrino_after_kf_E", &neutrino_after_kf_E, "neutrino_after_kf_E/F");

  float lepton_after_kf_pt, lepton_after_kf_eta, lepton_after_kf_phi, lepton_after_kf_et, lepton_after_kf_px, lepton_after_kf_py, lepton_after_kf_pz;
  tree->Branch("lepton_after_kf_pt", &lepton_after_kf_pt, "lepton_after_kf_pt/F");
  tree->Branch("lepton_after_kf_eta", &lepton_after_kf_eta, "lepton_after_kf_eta/F");
  tree->Branch("lepton_after_kf_phi", &lepton_after_kf_phi, "lepton_after_kf_phi/F");
  tree->Branch("lepton_after_kf_et", &lepton_after_kf_et, "lepton_after_kf_et/F");
  tree->Branch("lepton_after_kf_px", &lepton_after_kf_px, "lepton_after_kf_px/F");
  tree->Branch("lepton_after_kf_py", &lepton_after_kf_py, "lepton_after_kf_py/F");
  tree->Branch("lepton_after_kf_pz", &lepton_after_kf_pz, "lepton_after_kf_pz/F");
  float lepton_after_kf_E;
  tree->Branch("lepton_after_kf_E", &lepton_after_kf_E, "lepton_after_kf_E/F");

  float hadronic_top_mass_after_kf, leptonic_top_mass_after_kf, hadronic_w_mass_after_kf, leptonic_w_mass_after_kf;
  tree->Branch("hadronic_top_mass_after_kf", &hadronic_top_mass_after_kf, "hadronic_top_mass_after_kf/F");
  tree->Branch("leptonic_top_mass_after_kf", &leptonic_top_mass_after_kf, "leptonic_top_mass_after_kf/F");
  tree->Branch("hadronic_w_mass_after_kf", &hadronic_w_mass_after_kf, "hadronic_w_mass_after_kf/F");
  tree->Branch("leptonic_w_mass_after_kf", &leptonic_w_mass_after_kf, "leptonic_w_mass_after_kf/F");

  float kf_chisquare, kf_proba;
  tree->Branch("kf_chisquare", &kf_chisquare, "kf_chisquare/F");
  tree->Branch("kf_proba", &kf_proba, "kf_proba/F");

  tree->Branch("MC_mtt", &MC_mtt, "MC_mtt/F");

  // Gen informations
  float gen_lepton_pt, gen_lepton_eta, gen_lepton_phi, gen_lepton_et, gen_lepton_px, gen_lepton_py, gen_lepton_pz;
  tree->Branch("gen_lepton_pt", &gen_lepton_pt, "gen_lepton_pt/F");
  tree->Branch("gen_lepton_eta", &gen_lepton_eta, "gen_lepton_eta/F");
  tree->Branch("gen_lepton_phi", &gen_lepton_phi, "gen_lepton_phi/F");
  tree->Branch("gen_lepton_et", &gen_lepton_et, "gen_lepton_et/F");
  tree->Branch("gen_lepton_px", &gen_lepton_px, "gen_lepton_px/F");
  tree->Branch("gen_lepton_py", &gen_lepton_py, "gen_lepton_py/F");
  tree->Branch("gen_lepton_pz", &gen_lepton_pz, "gen_lepton_pz/F");
  float gen_lepton_E;
  tree->Branch("gen_lepton_E", &gen_lepton_E, "gen_lepton_E/F");

  float gen_firstJet_pt, gen_firstJet_eta, gen_firstJet_phi, gen_firstJet_et, gen_firstJet_px, gen_firstJet_py, gen_firstJet_pz;
  tree->Branch("gen_firstJet_pt", &gen_firstJet_pt, "gen_firstJet_pt/F");
  tree->Branch("gen_firstJet_eta", &gen_firstJet_eta, "gen_firstJet_eta/F");
  tree->Branch("gen_firstJet_phi", &gen_firstJet_phi, "gen_firstJet_phi/F");
  tree->Branch("gen_firstJet_et", &gen_firstJet_et, "gen_firstJet_et/F");
  tree->Branch("gen_firstJet_px", &gen_firstJet_px, "gen_firstJet_px/F");
  tree->Branch("gen_firstJet_py", &gen_firstJet_py, "gen_firstJet_py/F");
  tree->Branch("gen_firstJet_pz", &gen_firstJet_pz, "gen_firstJet_pz/F");
  float gen_firstJet_E;
  tree->Branch("gen_firstJet_E", &gen_firstJet_E, "gen_firstJet_E/F");

  float gen_secondJet_pt, gen_secondJet_eta, gen_secondJet_phi, gen_secondJet_et, gen_secondJet_px, gen_secondJet_py, gen_secondJet_pz;
  tree->Branch("gen_secondJet_pt", &gen_secondJet_pt, "gen_secondJet_pt/F");
  tree->Branch("gen_secondJet_eta", &gen_secondJet_eta, "gen_secondJet_eta/F");
  tree->Branch("gen_secondJet_phi", &gen_secondJet_phi, "gen_secondJet_phi/F");
  tree->Branch("gen_secondJet_et", &gen_secondJet_et, "gen_secondJet_et/F");
  tree->Branch("gen_secondJet_px", &gen_secondJet_px, "gen_secondJet_px/F");
  tree->Branch("gen_secondJet_py", &gen_secondJet_py, "gen_secondJet_py/F");
  tree->Branch("gen_secondJet_pz", &gen_secondJet_pz, "gen_secondJet_pz/F");
  float gen_secondJet_E;
  tree->Branch("gen_secondJet_E", &gen_secondJet_E, "gen_secondJet_E/F");

  float gen_hadronicB_pt, gen_hadronicB_eta, gen_hadronicB_phi, gen_hadronicB_et, gen_hadronicB_px, gen_hadronicB_py, gen_hadronicB_pz;
  tree->Branch("gen_hadronicB_pt", &gen_hadronicB_pt, "gen_hadronicB_pt/F");
  tree->Branch("gen_hadronicB_eta", &gen_hadronicB_eta, "gen_hadronicB_eta/F");
  tree->Branch("gen_hadronicB_phi", &gen_hadronicB_phi, "gen_hadronicB_phi/F");
  tree->Branch("gen_hadronicB_et", &gen_hadronicB_et, "gen_hadronicB_et/F");
  tree->Branch("gen_hadronicB_px", &gen_hadronicB_px, "gen_hadronicB_px/F");
  tree->Branch("gen_hadronicB_py", &gen_hadronicB_py, "gen_hadronicB_py/F");
  tree->Branch("gen_hadronicB_pz", &gen_hadronicB_pz, "gen_hadronicB_pz/F");
  float gen_hadronicB_E;
  tree->Branch("gen_hadronicB_E", &gen_hadronicB_E, "gen_hadronicB_E/F");

  float gen_leptonicB_pt, gen_leptonicB_eta, gen_leptonicB_phi, gen_leptonicB_et, gen_leptonicB_px, gen_leptonicB_py, gen_leptonicB_pz;
  tree->Branch("gen_leptonicB_pt", &gen_leptonicB_pt, "gen_leptonicB_pt/F");
  tree->Branch("gen_leptonicB_eta", &gen_leptonicB_eta, "gen_leptonicB_eta/F");
  tree->Branch("gen_leptonicB_phi", &gen_leptonicB_phi, "gen_leptonicB_phi/F");
  tree->Branch("gen_leptonicB_et", &gen_leptonicB_et, "gen_leptonicB_et/F");
  tree->Branch("gen_leptonicB_px", &gen_leptonicB_px, "gen_leptonicB_px/F");
  tree->Branch("gen_leptonicB_py", &gen_leptonicB_py, "gen_leptonicB_py/F");
  tree->Branch("gen_leptonicB_pz", &gen_leptonicB_pz, "gen_leptonicB_pz/F");
  float gen_leptonicB_E;
  tree->Branch("gen_leptonicB_E", &gen_leptonicB_E, "gen_leptonicB_E/F");

  float gen_neutrino_pt, gen_neutrino_eta, gen_neutrino_phi, gen_neutrino_et, gen_neutrino_px, gen_neutrino_py, gen_neutrino_pz;
  tree->Branch("gen_neutrino_pt", &gen_neutrino_pt, "gen_neutrino_pt/F");
  tree->Branch("gen_neutrino_eta", &gen_neutrino_eta, "gen_neutrino_eta/F");
  tree->Branch("gen_neutrino_phi", &gen_neutrino_phi, "gen_neutrino_phi/F");
  tree->Branch("gen_neutrino_et", &gen_neutrino_et, "gen_neutrino_et/F");
  tree->Branch("gen_neutrino_px", &gen_neutrino_px, "gen_neutrino_px/F");
  tree->Branch("gen_neutrino_py", &gen_neutrino_py, "gen_neutrino_py/F");
  tree->Branch("gen_neutrino_pz", &gen_neutrino_pz, "gen_neutrino_pz/F");
  float gen_neutrino_E;
  tree->Branch("gen_neutrino_E", &gen_neutrino_E, "gen_neutrino_E/F");

  TProfile* mtt_reco_vs_mtt = new TProfile("mtt_reco_vs_mtt_profile", ";Generated m_{t#bar{t}};m_{t#bar{t}} reco", nBins, bins);
  TProfile* mtt_reco_kf_vs_mtt = new TProfile("mtt_reco_kf_vs_mtt_profile", ";Generated m_{t#bar{t}};m_{t#bar{t}} reco after KinFit", nBins, bins);

  GaussianProfile* h_mtt_reco_vs_mtt = new GaussianProfile("mtt_reco_vs_mtt", nBins, bins);
  h_mtt_reco_vs_mtt->setPrefix("mtt");

  //std::vector<TH1F*> h_mtt_reco_vs_mtt = buildVector<TH1F>("mtt_reco_vs_mtt", 0, 0, 0);
  GaussianProfile* h_hadronic_top_reco_vs_mtt = new GaussianProfile("hadronic_top_reco_vs_mtt", nBins, bins, 50, 140, 250);
  GaussianProfile* h_leptonic_top_reco_vs_mtt = new GaussianProfile("leptonic_top_reco_vs_mtt", nBins, bins, 50, 140, 250);
  GaussianProfile* h_hadronic_w_reco_vs_mtt = new GaussianProfile("hadronic_w_reco_vs_mtt", nBins, bins, 50, 40, 120);
  GaussianProfile* h_leptonic_w_reco_vs_mtt = new GaussianProfile("leptonic_w_reco_vs_mtt", nBins, bins, 50, 40, 120);

  GaussianProfile* h_mtt_reco_after_kf_vs_mtt = new GaussianProfile("mtt_reco_after_kf_vs_mtt", nBins, bins);

  GaussianProfile* h_hadronic_top_reco_after_kf_vs_mtt = new GaussianProfile("hadronic_top_reco_after_kf_vs_mtt", nBins, bins, 50, 140, 250);
  GaussianProfile* h_leptonic_top_reco_after_kf_vs_mtt = new GaussianProfile("leptonic_top_reco_after_kf_vs_mtt", nBins, bins, 50, 140, 250);
  GaussianProfile* h_hadronic_w_reco_after_kf_vs_mtt = new GaussianProfile("hadronic_w_reco_after_kf_vs_mtt", nBins, bins, 50, 40, 120);
  GaussianProfile* h_leptonic_w_reco_after_kf_vs_mtt = new GaussianProfile("leptonic_w_reco_after_kf_vs_mtt", nBins, bins, 50, 40, 120);

  GaussianProfile* h_delta_leptonicB_Et_vs_mtt = new GaussianProfile("delta_leptonic_B_Et_vs_mtt", nBins, bins, 100, -200, 200, false);
  GaussianProfile* h_delta_leptonicB_Eta_vs_mtt = new GaussianProfile("delta_leptonic_B_Eta_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);
  GaussianProfile* h_delta_leptonicB_Phi_vs_mtt = new GaussianProfile("delta_leptonic_B_Phi_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);

  GaussianProfile* h_delta_hadronicB_Et_vs_mtt = new GaussianProfile("delta_hadronic_B_Et_vs_mtt", nBins, bins, 100, -200, 200, false);
  GaussianProfile* h_delta_hadronicB_Eta_vs_mtt = new GaussianProfile("delta_hadronic_B_Eta_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);
  GaussianProfile* h_delta_hadronicB_Phi_vs_mtt = new GaussianProfile("delta_hadronic_B_Phi_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);

  GaussianProfile* h_delta_firstJet_Et_vs_mtt = new GaussianProfile("delta_firstJet_Et_vs_mtt", nBins, bins, 100, -200, 200, false);
  GaussianProfile* h_delta_firstJet_Eta_vs_mtt = new GaussianProfile("delta_firstJet_Eta_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);
  GaussianProfile* h_delta_firstJet_Phi_vs_mtt = new GaussianProfile("delta_firstJet_Phi_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);

  GaussianProfile* h_delta_secondJet_Et_vs_mtt = new GaussianProfile("delta_secondJet_Et_vs_mtt", nBins, bins, 100, -200, 200, false);
  GaussianProfile* h_delta_secondJet_Eta_vs_mtt = new GaussianProfile("delta_secondJet_Eta_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);
  GaussianProfile* h_delta_secondJet_Phi_vs_mtt = new GaussianProfile("delta_secondJet_Phi_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);

  GaussianProfile* h_delta_lepton_Et_vs_mtt = new GaussianProfile("delta_lepton_Et_vs_mtt", nBins, bins, 100, -200, 200, false);
  GaussianProfile* h_delta_lepton_Eta_vs_mtt = new GaussianProfile("delta_lepton_Eta_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);
  GaussianProfile* h_delta_lepton_Phi_vs_mtt = new GaussianProfile("delta_lepton_Phi_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);

  GaussianProfile* h_delta_neutrino_Et_vs_mtt = new GaussianProfile("delta_neutrino_Et_vs_mtt", nBins, bins, 100, -200, 200, false);
  GaussianProfile* h_delta_neutrino_Eta_vs_mtt = new GaussianProfile("delta_neutrino_Eta_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);
  GaussianProfile* h_delta_neutrino_Phi_vs_mtt = new GaussianProfile("delta_neutrino_Phi_vs_mtt", nBins, bins, 100, -3.14, 3.14, false);

  TChain* MC = NULL, *event = NULL, *jets = NULL, *MET = NULL, *muons = NULL, *electrons = NULL, *mtt = NULL;
  loadChain(inputFiles, MC, event, jets, MET, muons, electrons, mtt);

  MC->SetBranchStatus("*", 0);
  event->SetBranchStatus("*", 0);
  jets->SetBranchStatus("*", 0);
  MET->SetBranchStatus("*", 0);
  muons->SetBranchStatus("*", 0);
  electrons->SetBranchStatus("*", 0);
  mtt->SetBranchStatus("*", 0);

  TClonesArray *m_MC_lorentzvector;
  int   m_n_MCs;
  m_MC_lorentzvector = new TClonesArray("TLorentzVector");

  SetBranchAddress(MC, "n_MCs",  &m_n_MCs);
  SetBranchAddress(MC, "MC_4vector", &m_MC_lorentzvector);

  // Jets
  SetBranchAddress(jets, "n_jets", &n_jets);
  SetBranchAddress(jets, "jet_4vector", &m_jet_lorentzvector);
  SetBranchAddress(jets, "jet_mcParticleIndex", &m_jet_MCIndex);

  // MTT
  SetBranchAddress(mtt, "MC_leptonicBIndex", &genLeptonicBIndex);
  SetBranchAddress(mtt, "MC_hadronicBIndex", &genHadronicBIndex);
  SetBranchAddress(mtt, "MC_hadronicFirstJetIndex", &genFirstJetIndex);
  SetBranchAddress(mtt, "MC_hadronicSecondJetIndex", &genSecondJetIndex);
  SetBranchAddress(mtt, "MC_leptonIndex", &genLeptonIndex);
  SetBranchAddress(mtt, "MC_neutrinoIndex", &genNeutrinoIndex);

  /*
  int selectedLeptonicBIndex, selectedHadronicBIndex, selectedFirstJetIndex, selectedSecondJetIndex;
  SetBranchAddress(mtt, "selectedLeptonicBIndex", &selectedLeptonicBIndex);
  SetBranchAddress(mtt, "selectedHadronicBIndex", &selectedHadronicBIndex);
  SetBranchAddress(mtt, "selectedHadronicFirstJetIndex", &selectedFirstJetIndex);
  SetBranchAddress(mtt, "selectedHadronicSecondJetIndex", &selectedSecondJetIndex);
  */

  int numComb, channel;
  SetBranchAddress(mtt, "MC_channel", &channel);
  //SetBranchAddress(mtt, "isSel", &isSel);
  SetBranchAddress(mtt, "numComb", &numComb);
  SetBranchAddress(mtt, "MC_mtt", &MC_mtt);

    // Electrons
  static const int 	m_electrons_MAX  = 100;
  TClonesArray* m_ele_lorentzvector = new TClonesArray("TLorentzVector");
  int n_electrons;
  int m_ele_MCIndex[m_electrons_MAX];
  int ele_charge[m_electrons_MAX];
  
  SetBranchAddress(electrons, "n_electrons", &n_electrons);
  SetBranchAddress(electrons, "electron_4vector", &m_ele_lorentzvector);
  SetBranchAddress(electrons, "electron_mcParticleIndex", &m_ele_MCIndex);
  SetBranchAddress(electrons, "electron_charge",          &ele_charge);

  // Muons
  static const int 	m_muons_MAX  = 100;
  TClonesArray* m_muo_lorentzvector = new TClonesArray("TLorentzVector");
  int n_muons;
  int m_muo_MCIndex[m_muons_MAX];
  int muon_charge[m_muons_MAX];

  SetBranchAddress(muons, "n_muons",  &n_muons);
  SetBranchAddress(muons, "muon_4vector", &m_muo_lorentzvector);
  SetBranchAddress(muons, "muon_mcParticleIndex", &m_muo_MCIndex);
  SetBranchAddress(muons, "muon_charge",          &muon_charge);

  // MET
  TClonesArray* m_met_lorentzvector = new TClonesArray("TLorentzVector");
  SetBranchAddress(MET, "met_4vector", &m_met_lorentzvector);

  uint64_t possibleMatchableEntries = 0;

  uint64_t totalFittedEntries = 0;
  uint64_t totalFittedConvergedEntries = 0;

  uint64_t entries = MC->GetEntries();
  //uint64_t entries = 10000;
  if (ENTRIES > 0 && ENTRIES <= entries)
    entries = ENTRIES;

  for (uint64_t entry = 0; entry < entries; entry++) {
    MC->GetEntry(entry);
    event->GetEntry(entry);
    jets->GetEntry(entry);
    MET->GetEntry(entry);
    muons->GetEntry(entry);
    electrons->GetEntry(entry);
    mtt->GetEntry(entry);

    if (((entry) % 100000) == 0) {
      std::cout << "Processing entry " << entry + 1 << " out of " << entries << std::endl;
    }
    /*if (((entry) % 1) == 0) {
      std::cout << "Processing entry " << entry + 1 << " out of " << entries << std::endl;
    }*/

    if (genLeptonicBIndex == -1 || genHadronicBIndex == -1 || genFirstJetIndex == -1 || genHadronicBIndex == -1)
      continue;

    if (channel != 1 && channel != 2)
      continue;
    
    if (numComb) {
      possibleMatchableEntries++;
    }

    // Ensure gen particles are matched to one reco jets
    if (!isMatched(genLeptonicBIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets) ||
        !isMatched(genHadronicBIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets) ||
        !isMatched(genFirstJetIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets)  ||
        !isMatched(genSecondJetIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets))
      continue;

    // MET
    TLorentzVector* neutrinoP4 = static_cast<TLorentzVector*>((*m_met_lorentzvector)[0]);
    if (neutrinoP4->Pt() < 20)
      continue;

    // Lepton
    TLorentzVector* leptonP4 = NULL;
    int leptonCharge = 0;
    float ptLeptonCut = 0;
    if (channel == 2) {
      if (n_muons != 1)
        continue;

      leptonP4 = static_cast<TLorentzVector*>((*m_muo_lorentzvector)[0]);
      leptonCharge = muon_charge[0];
      ptLeptonCut = 25;
    } else {
      if (n_electrons != 1)
        continue;

      leptonP4 = static_cast<TLorentzVector*>((*m_ele_lorentzvector)[0]);
      leptonCharge = ele_charge[0];
      ptLeptonCut = 30;
    }

    if (leptonP4->Pt() < ptLeptonCut || fabs(leptonP4->Eta()) > 2.1)
      continue;

    TLorentzVector* firstJetP4 = NULL, *secondJetP4 = NULL, *leptonicBP4 = NULL, *hadronicBP4 = NULL;

    for (int i = 0; i < n_jets; i++) {
      int index = m_jet_MCIndex[i];
      if (index == genLeptonicBIndex)
        leptonicBP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[i]);
      else if (index == genHadronicBIndex)
        hadronicBP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[i]);
      else if (index == genFirstJetIndex)
        firstJetP4  = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[i]);
      else if (index == genSecondJetIndex)
        secondJetP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[i]);
    }

    std::vector<double> pts = {firstJetP4->Pt(), secondJetP4->Pt(), hadronicBP4->Pt(), leptonicBP4->Pt()};
    std::sort(pts.begin(), pts.end());

    if (pts[3] < 70 || pts[2] < 50 || pts[1] < 30 || pts[0] < 30)
      continue;

    if (fabs(firstJetP4->Eta()) > 2.4 || fabs(secondJetP4->Eta()) > 2.4 || fabs(hadronicBP4->Eta()) > 2.4 || fabs(leptonicBP4->Eta()) > 2.4)
      continue;

    if (! kinFitter.PzNeutrino(*leptonP4, *neutrinoP4, *leptonicBP4))
      continue;

    TLorentzVector* gen_leptonicBP4 = static_cast<TLorentzVector*>((*m_MC_lorentzvector)[genLeptonicBIndex]);
    TLorentzVector* gen_hadronicBP4 = static_cast<TLorentzVector*>((*m_MC_lorentzvector)[genHadronicBIndex]);
    TLorentzVector* gen_firstJetP4 = static_cast<TLorentzVector*>((*m_MC_lorentzvector)[genFirstJetIndex]);
    TLorentzVector* gen_secondJetP4 = static_cast<TLorentzVector*>((*m_MC_lorentzvector)[genSecondJetIndex]);
    TLorentzVector* gen_leptonP4 = static_cast<TLorentzVector*>((*m_MC_lorentzvector)[genLeptonIndex]);
    TLorentzVector* gen_neutrinoP4 = static_cast<TLorentzVector*>((*m_MC_lorentzvector)[genNeutrinoIndex]);

    mtt_reco = (*firstJetP4 + *secondJetP4 + *hadronicBP4 + *leptonicBP4 + *leptonP4 + *neutrinoP4).M();
    hadronic_top_mass = (*firstJetP4 + *secondJetP4 + *hadronicBP4).M();
    leptonic_top_mass = (*leptonicBP4 + *leptonP4 + *neutrinoP4).M();
    hadronic_w_mass = (*firstJetP4 + *secondJetP4).M();
    leptonic_w_mass = (*leptonP4 + *neutrinoP4).M();

    mtt_reco_vs_mtt->Fill(mtt_reco, MC_mtt);

    h_mtt_reco_vs_mtt->fill(MC_mtt, mtt_reco);
    h_hadronic_top_reco_vs_mtt->fill(MC_mtt, hadronic_top_mass);
    h_leptonic_top_reco_vs_mtt->fill(MC_mtt, leptonic_top_mass);
    h_hadronic_w_reco_vs_mtt->fill(MC_mtt, hadronic_w_mass);
    h_leptonic_w_reco_vs_mtt->fill(MC_mtt, leptonic_w_mass);

    {
      // Fill tree variables
      firstJet_pt = firstJetP4->Pt();
      firstJet_eta = firstJetP4->Eta();
      firstJet_phi = firstJetP4->Phi();
      firstJet_et  = firstJetP4->Et();
      firstJet_px  = firstJetP4->Px();
      firstJet_py  = firstJetP4->Py();
      firstJet_pz  = firstJetP4->Pz();
      firstJet_E   = firstJetP4->E();

      secondJet_pt = secondJetP4->Pt();
      secondJet_eta = secondJetP4->Eta();
      secondJet_phi = secondJetP4->Phi();
      secondJet_et  = secondJetP4->Et();
      secondJet_px  = secondJetP4->Px();
      secondJet_py  = secondJetP4->Py();
      secondJet_pz  = secondJetP4->Pz();
      secondJet_E   = secondJetP4->E();

      hadronicB_pt = hadronicBP4->Pt();
      hadronicB_eta = hadronicBP4->Eta();
      hadronicB_phi = hadronicBP4->Phi();
      hadronicB_et  = hadronicBP4->Et();
      hadronicB_px  = hadronicBP4->Px();
      hadronicB_py  = hadronicBP4->Py();
      hadronicB_pz  = hadronicBP4->Pz();
      hadronicB_E   = hadronicBP4->E();

      leptonicB_pt = leptonicBP4->Pt();
      leptonicB_eta = leptonicBP4->Eta();
      leptonicB_phi = leptonicBP4->Phi();
      leptonicB_et  = leptonicBP4->Et();
      leptonicB_px  = leptonicBP4->Px();
      leptonicB_py  = leptonicBP4->Py();
      leptonicB_pz  = leptonicBP4->Pz();
      leptonicB_E   = leptonicBP4->E();

      neutrino_pt = neutrinoP4->Pt();
      neutrino_eta = neutrinoP4->Eta();
      neutrino_phi = neutrinoP4->Phi();
      neutrino_et  = neutrinoP4->Et();
      neutrino_px  = neutrinoP4->Px();
      neutrino_py  = neutrinoP4->Py();
      neutrino_pz  = neutrinoP4->Pz();
      neutrino_E   = neutrinoP4->E();

      lepton_pt = leptonP4->Pt();
      lepton_eta = leptonP4->Eta();
      lepton_phi = leptonP4->Phi();
      lepton_et  = leptonP4->Et();
      lepton_px  = leptonP4->Px();
      lepton_py  = leptonP4->Py();
      lepton_pz  = leptonP4->Pz();
      lepton_E   = leptonP4->E();

      gen_firstJet_pt = gen_firstJetP4->Pt();
      gen_firstJet_eta = gen_firstJetP4->Eta();
      gen_firstJet_phi = gen_firstJetP4->Phi();
      gen_firstJet_et  = gen_firstJetP4->Et();
      gen_firstJet_px  = gen_firstJetP4->Px();
      gen_firstJet_py  = gen_firstJetP4->Py();
      gen_firstJet_pz  = gen_firstJetP4->Pz();
      gen_firstJet_E   = gen_firstJetP4->E();

      gen_secondJet_pt = gen_secondJetP4->Pt();
      gen_secondJet_eta = gen_secondJetP4->Eta();
      gen_secondJet_phi = gen_secondJetP4->Phi();
      gen_secondJet_et  = gen_secondJetP4->Et();
      gen_secondJet_px  = gen_secondJetP4->Px();
      gen_secondJet_py  = gen_secondJetP4->Py();
      gen_secondJet_pz  = gen_secondJetP4->Pz();
      gen_secondJet_E   = gen_secondJetP4->E();

      gen_hadronicB_pt = gen_hadronicBP4->Pt();
      gen_hadronicB_eta = gen_hadronicBP4->Eta();
      gen_hadronicB_phi = gen_hadronicBP4->Phi();
      gen_hadronicB_et  = gen_hadronicBP4->Et();
      gen_hadronicB_px  = gen_hadronicBP4->Px();
      gen_hadronicB_py  = gen_hadronicBP4->Py();
      gen_hadronicB_pz  = gen_hadronicBP4->Pz();
      gen_hadronicB_E   = gen_hadronicBP4->E();

      gen_leptonicB_pt = gen_leptonicBP4->Pt();
      gen_leptonicB_eta = gen_leptonicBP4->Eta();
      gen_leptonicB_phi = gen_leptonicBP4->Phi();
      gen_leptonicB_et  = gen_leptonicBP4->Et();
      gen_leptonicB_px  = gen_leptonicBP4->Px();
      gen_leptonicB_py  = gen_leptonicBP4->Py();
      gen_leptonicB_pz  = gen_leptonicBP4->Pz();
      gen_leptonicB_E   = gen_leptonicBP4->E();

      gen_neutrino_pt = gen_neutrinoP4->Pt();
      gen_neutrino_eta = gen_neutrinoP4->Eta();
      gen_neutrino_phi = gen_neutrinoP4->Phi();
      gen_neutrino_et  = gen_neutrinoP4->Et();
      gen_neutrino_px  = gen_neutrinoP4->Px();
      gen_neutrino_py  = gen_neutrinoP4->Py();
      gen_neutrino_pz  = gen_neutrinoP4->Pz();
      gen_neutrino_E   = gen_neutrinoP4->E();

      gen_lepton_pt = gen_leptonP4->Pt();
      gen_lepton_eta = gen_leptonP4->Eta();
      gen_lepton_phi = gen_leptonP4->Phi();
      gen_lepton_et  = gen_leptonP4->Et();
      gen_lepton_px  = gen_leptonP4->Px();
      gen_lepton_py  = gen_leptonP4->Py();
      gen_lepton_pz  = gen_leptonP4->Pz();
      gen_lepton_E   = gen_leptonP4->E();
    }
      

    totalFittedEntries++;

    //std::cout << "Before fit" << std::endl;
    if (fitter.fit(*firstJetP4, *secondJetP4, *hadronicBP4, *leptonicBP4, *leptonP4, *neutrinoP4, leptonCharge, (channel == 1) ? CovarianceMatrix::kElectron : CovarianceMatrix::kMuon) != 0) {
      //std::cout << "Fit failed!" << std::endl;
      continue;
    }

    totalFittedConvergedEntries++;

    //std::cout << "After fit" << std::endl;

    //std::cout << "Chi2: " << fitter.fitS() << std::endl;
    kf_chisquare = fitter.fitS();
    kf_proba     = fitter.fitProb();

    /*firstJetP4 = kinFit->GetFittedFirstLightJet();
    secondJetP4 = kinFit->GetFittedSecondLightJet();
    hadronicBP4 = kinFit->GetFittedHadronicBJet();
    leptonicBP4 = kinFit->GetFittedLeptonicBJet();
    leptonP4 = kinFit->GetFittedLepton();
    neutrinoP4 = kinFit->GetFittedNeutrino();*/

    const reco::Candidate::PolarLorentzVector firstJetP4_after_kf = fitter.fittedHadP().polarP4();
    const reco::Candidate::PolarLorentzVector secondJetP4_after_kf = fitter.fittedHadQ().polarP4();

    const reco::Candidate::PolarLorentzVector hadronicBP4_after_kf = fitter.fittedHadB().polarP4();
    const reco::Candidate::PolarLorentzVector leptonicBP4_after_kf = fitter.fittedLepB().polarP4();

    const reco::Candidate::PolarLorentzVector leptonP4_after_kf = fitter.fittedLepton().polarP4();
    const reco::Candidate::PolarLorentzVector neutrinoP4_after_kf = fitter.fittedNeutrino().polarP4();

    mtt_reco_after_kf = (firstJetP4_after_kf + secondJetP4_after_kf + hadronicBP4_after_kf + leptonicBP4_after_kf + leptonP4_after_kf + neutrinoP4_after_kf).M();
    hadronic_top_mass_after_kf = (firstJetP4_after_kf + secondJetP4_after_kf + hadronicBP4_after_kf).M();
    leptonic_top_mass_after_kf = (leptonicBP4_after_kf + leptonP4_after_kf + neutrinoP4_after_kf).M();
    hadronic_w_mass_after_kf = (firstJetP4_after_kf + secondJetP4_after_kf).M();
    leptonic_w_mass_after_kf = (leptonP4_after_kf + neutrinoP4_after_kf).M();

    mtt_reco_kf_vs_mtt->Fill(mtt_reco_after_kf, MC_mtt);
    h_mtt_reco_after_kf_vs_mtt->fill(MC_mtt, mtt_reco_after_kf);
    h_hadronic_top_reco_after_kf_vs_mtt->fill(MC_mtt, hadronic_top_mass_after_kf);
    h_leptonic_top_reco_after_kf_vs_mtt->fill(MC_mtt, leptonic_top_mass_after_kf);
    h_hadronic_w_reco_after_kf_vs_mtt->fill(MC_mtt, hadronic_w_mass_after_kf);
    h_leptonic_w_reco_after_kf_vs_mtt->fill(MC_mtt, leptonic_w_mass_after_kf);

    h_delta_leptonicB_Et_vs_mtt->fill(MC_mtt, leptonicBP4_after_kf.Et() - leptonicBP4->Et());
    h_delta_leptonicB_Eta_vs_mtt->fill(MC_mtt, leptonicBP4_after_kf.Eta() - leptonicBP4->Eta());
    h_delta_leptonicB_Phi_vs_mtt->fill(MC_mtt, leptonicBP4_after_kf.Phi() - leptonicBP4->Phi());

    h_delta_hadronicB_Et_vs_mtt->fill(MC_mtt, hadronicBP4_after_kf.Et() - hadronicBP4->Et());
    h_delta_hadronicB_Eta_vs_mtt->fill(MC_mtt, hadronicBP4_after_kf.Eta() - hadronicBP4->Eta());
    h_delta_hadronicB_Phi_vs_mtt->fill(MC_mtt, hadronicBP4_after_kf.Phi() - hadronicBP4->Phi());

    h_delta_firstJet_Et_vs_mtt->fill(MC_mtt, firstJetP4_after_kf.Et() - firstJetP4->Et());
    h_delta_firstJet_Eta_vs_mtt->fill(MC_mtt, firstJetP4_after_kf.Eta() - firstJetP4->Eta());
    h_delta_firstJet_Phi_vs_mtt->fill(MC_mtt, firstJetP4_after_kf.Phi() - firstJetP4->Phi());

    h_delta_secondJet_Et_vs_mtt->fill(MC_mtt, secondJetP4_after_kf.Et() - secondJetP4->Et());
    h_delta_secondJet_Eta_vs_mtt->fill(MC_mtt, secondJetP4_after_kf.Eta() - secondJetP4->Eta());
    h_delta_secondJet_Phi_vs_mtt->fill(MC_mtt, secondJetP4_after_kf.Phi() - secondJetP4->Phi());

    h_delta_lepton_Et_vs_mtt->fill(MC_mtt, leptonP4_after_kf.Et() - leptonP4->Et());
    h_delta_lepton_Eta_vs_mtt->fill(MC_mtt, leptonP4_after_kf.Eta() - leptonP4->Eta());
    h_delta_lepton_Phi_vs_mtt->fill(MC_mtt, leptonP4_after_kf.Phi() - leptonP4->Phi());

    h_delta_neutrino_Et_vs_mtt->fill(MC_mtt, neutrinoP4_after_kf.Et() - neutrinoP4->Et());
    h_delta_neutrino_Eta_vs_mtt->fill(MC_mtt, neutrinoP4_after_kf.Eta() - neutrinoP4->Eta());
    h_delta_neutrino_Phi_vs_mtt->fill(MC_mtt, neutrinoP4_after_kf.Phi() - neutrinoP4->Phi());

    //std::cout << hadronic_top_mass_after_kf << std::endl;

    {
      // Fill tree variables
      firstJet_after_kf_pt = firstJetP4_after_kf.Pt();
      firstJet_after_kf_eta = firstJetP4_after_kf.Eta();
      firstJet_after_kf_phi = firstJetP4_after_kf.Phi();
      firstJet_after_kf_et  = firstJetP4_after_kf.Et();
      firstJet_after_kf_px  = firstJetP4_after_kf.Px();
      firstJet_after_kf_py  = firstJetP4_after_kf.Py();
      firstJet_after_kf_pz  = firstJetP4_after_kf.Pz();
      firstJet_after_kf_E   = firstJetP4_after_kf.E();

      secondJet_after_kf_pt = secondJetP4_after_kf.Pt();
      secondJet_after_kf_eta = secondJetP4_after_kf.Eta();
      secondJet_after_kf_phi = secondJetP4_after_kf.Phi();
      secondJet_after_kf_et  = secondJetP4_after_kf.Et();
      secondJet_after_kf_px  = secondJetP4_after_kf.Px();
      secondJet_after_kf_py  = secondJetP4_after_kf.Py();
      secondJet_after_kf_pz  = secondJetP4_after_kf.Pz();
      secondJet_after_kf_E   = secondJetP4_after_kf.E();

      hadronicB_after_kf_pt = hadronicBP4_after_kf.Pt();
      hadronicB_after_kf_eta = hadronicBP4_after_kf.Eta();
      hadronicB_after_kf_phi = hadronicBP4_after_kf.Phi();
      hadronicB_after_kf_et  = hadronicBP4_after_kf.Et();
      hadronicB_after_kf_px  = hadronicBP4_after_kf.Px();
      hadronicB_after_kf_py  = hadronicBP4_after_kf.Py();
      hadronicB_after_kf_pz  = hadronicBP4_after_kf.Pz();
      hadronicB_after_kf_E   = hadronicBP4_after_kf.E();

      leptonicB_after_kf_pt = leptonicBP4_after_kf.Pt();
      leptonicB_after_kf_eta = leptonicBP4_after_kf.Eta();
      leptonicB_after_kf_phi = leptonicBP4_after_kf.Phi();
      leptonicB_after_kf_et  = leptonicBP4_after_kf.Et();
      leptonicB_after_kf_px  = leptonicBP4_after_kf.Px();
      leptonicB_after_kf_py  = leptonicBP4_after_kf.Py();
      leptonicB_after_kf_pz  = leptonicBP4_after_kf.Pz();
      leptonicB_after_kf_E   = leptonicBP4_after_kf.E();

      neutrino_after_kf_pt = neutrinoP4_after_kf.Pt();
      neutrino_after_kf_eta = neutrinoP4_after_kf.Eta();
      neutrino_after_kf_phi = neutrinoP4_after_kf.Phi();
      neutrino_after_kf_et  = neutrinoP4_after_kf.Et();
      neutrino_after_kf_px  = neutrinoP4_after_kf.Px();
      neutrino_after_kf_py  = neutrinoP4_after_kf.Py();
      neutrino_after_kf_pz  = neutrinoP4_after_kf.Pz();
      neutrino_after_kf_E   = neutrinoP4_after_kf.E();

      lepton_after_kf_pt = leptonP4_after_kf.Pt();
      lepton_after_kf_eta = leptonP4_after_kf.Eta();
      lepton_after_kf_phi = leptonP4_after_kf.Phi();
      lepton_after_kf_et  = leptonP4_after_kf.Et();
      lepton_after_kf_px  = leptonP4_after_kf.Px();
      lepton_after_kf_py  = leptonP4_after_kf.Py();
      lepton_after_kf_pz  = leptonP4_after_kf.Pz();
      lepton_after_kf_E   = leptonP4_after_kf.E();
    }

    tree->Fill();
  }

  std::cout << "KinFit efficiency: " << (double) totalFittedConvergedEntries / totalFittedEntries * 100 << " % (" << totalFittedConvergedEntries << " / " << totalFittedEntries << ")" << std::endl;

  delete MC; delete event; delete MET; delete jets; delete muons; delete electrons;

  if (outputFile.length() > 0) {
    TFile* f = TFile::Open(outputFile.c_str(), "recreate");
    mtt_reco_vs_mtt->Write();
    mtt_reco_kf_vs_mtt->Write();

    h_mtt_reco_vs_mtt->write(f);
    h_hadronic_top_reco_vs_mtt->write(f);
    h_leptonic_top_reco_vs_mtt->write(f);
    h_hadronic_w_reco_vs_mtt->write(f);
    h_leptonic_w_reco_vs_mtt->write(f);

    h_mtt_reco_after_kf_vs_mtt->write(f);
    h_hadronic_top_reco_after_kf_vs_mtt->write(f);
    h_leptonic_top_reco_after_kf_vs_mtt->write(f);
    h_hadronic_w_reco_after_kf_vs_mtt->write(f);
    h_leptonic_w_reco_after_kf_vs_mtt->write(f);

    h_delta_leptonicB_Et_vs_mtt->write(f);
    h_delta_leptonicB_Eta_vs_mtt->write(f);
    h_delta_leptonicB_Phi_vs_mtt->write(f);

    h_delta_hadronicB_Et_vs_mtt->write(f);
    h_delta_hadronicB_Eta_vs_mtt->write(f);
    h_delta_hadronicB_Phi_vs_mtt->write(f);

    h_delta_firstJet_Et_vs_mtt->write(f);
    h_delta_firstJet_Eta_vs_mtt->write(f);
    h_delta_firstJet_Phi_vs_mtt->write(f);

    h_delta_secondJet_Et_vs_mtt->write(f);
    h_delta_secondJet_Eta_vs_mtt->write(f);
    h_delta_secondJet_Phi_vs_mtt->write(f);

    h_delta_lepton_Et_vs_mtt->write(f);
    h_delta_lepton_Eta_vs_mtt->write(f);
    h_delta_lepton_Phi_vs_mtt->write(f);

    h_delta_neutrino_Et_vs_mtt->write(f);
    h_delta_neutrino_Eta_vs_mtt->write(f);
    h_delta_neutrino_Phi_vs_mtt->write(f);

    tree->Write();

    f->Close();
    delete f;
  }

  delete h_mtt_reco_vs_mtt;
  delete h_hadronic_top_reco_vs_mtt;
  delete h_leptonic_top_reco_vs_mtt;
  delete h_hadronic_w_reco_vs_mtt;
  delete h_leptonic_w_reco_vs_mtt;
  delete h_mtt_reco_after_kf_vs_mtt;
  delete h_hadronic_top_reco_after_kf_vs_mtt;
  delete h_leptonic_top_reco_after_kf_vs_mtt;
  delete h_hadronic_w_reco_after_kf_vs_mtt;
  delete h_leptonic_w_reco_after_kf_vs_mtt;

  delete h_delta_leptonicB_Et_vs_mtt;
  delete h_delta_leptonicB_Eta_vs_mtt;
  delete h_delta_leptonicB_Phi_vs_mtt;

  delete h_delta_hadronicB_Et_vs_mtt;
  delete h_delta_hadronicB_Eta_vs_mtt;
  delete h_delta_hadronicB_Phi_vs_mtt;

  delete tree;

  //delete kinFit;
}

void loadInputFiles(const std::string& filename, std::vector<std::string>& files) {

  ifstream ifs(filename.c_str());
  std::string line;

  while (getline(ifs, line))
    files.push_back(line);

  ifs.close();
}

int main(int argc, char** argv)
{
  //gErrorIgnoreLevel = kFatal;

  //edm::service::MessageServicePresence my_message_service;

  try {
    TCLAP::CmdLine cmd("compute various distribution to optimize chi square", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> configFileArg("c", "cfg", "Python config file", true, "", "string", cmd);

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", false, "", "string", cmd);

    TCLAP::ValueArg<uint64_t> entriesArg("n", "entries", "Number of entries to process", false, 0, "int", cmd);

    /*TCLAP::ValueArg<std::string> typeArg("", "type", "current inputfile type (semie or semimu)", false, "", "string", cmd);
      TCLAP::ValueArg<std::string> pileupArg("", "pileup", "PU profile used for MC production", false, "S7", "string", cmd);*/

    TCLAP::SwitchArg csvModeArg("", "csv", "csv mode", cmd);

    cmd.parse(argc, argv);

    // Parse python config file
    boost::shared_ptr<edm::ParameterSet> rootPS = edm::readPSetsFrom(configFileArg.getValue());
    if (! rootPS->existsAs<edm::ParameterSet>("process") ){
      std::cout << " ERROR: ParametersSet 'process' is missing in your configuration file" << std::endl;
      return 0;
    }

    const edm::ParameterSet& psProcess = rootPS->getParameter<edm::ParameterSet>("process");
    const edm::ParameterSet& kinFitPS = psProcess.getParameter<edm::ParameterSet>("kinfit");

    //std::cout << kinFitPS.dump() << std::endl;

    CSV_MODE = csvModeArg.getValue();
    ENTRIES = entriesArg.getValue();

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    process(inputFiles, outputFileArg.getValue(), kinFitPS);    

  } catch (TCLAP::ArgException& e) {
    std::cout << e.what() << std::endl;
  }

}
