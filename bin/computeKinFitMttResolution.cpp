/**
 * Based on assiciable events, run the KinFit on all jets combinaisons, and
 * find the number of time the right solution (good jets placement) is obtained
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

int genLeptonicBIndex, genHadronicBIndex, genFirstJetIndex, genSecondJetIndex;

const int nBins = 15;
const double bins[] = {340, 360, 380, 400, 420, 460, 500, 550, 600, 650, 750, 850, 950, 1050, 1200, 1400, 1600};

bool CSV_MODE = true;
bool SEMI_MU  = false;
uint64_t ENTRIES = 0;

double binomialError(double efficiency, int n) {
  return sqrt((efficiency * (1 - efficiency)) / (n));
}

void printEff(double efficiency, int n) {
  if (CSV_MODE) {
    std::cout << std::left << std::setw(15) << efficiency << "\t" << binomialError(efficiency, n) << std::endl;
  } else {
    std::cout << efficiency * 100 << " +/- " << binomialError(efficiency, n) * 100 << std::endl;
  }
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

  float kf_chisquare, kf_proba;
  tree->Branch("kf_chisquare", &kf_chisquare, "kf_chisquare/F");
  tree->Branch("kf_proba", &kf_proba, "kf_proba/F");

  tree->Branch("MC_mtt", &MC_mtt, "MC_mtt/F");

  bool goodCombinaison = false;
  tree->Branch("good_combinaison", &goodCombinaison, "good_combinaison/O");

  bool wellPlacedAfterKF = false;
  tree->Branch("combinaison_well_placed_after_kf", &wellPlacedAfterKF, "combinaison_well_placed_after_kf/O");

  bool associable = false;
  tree->Branch("associable", &associable, "associable/O");

  bool kf_converged = false;
  tree->Branch("kf_converged", &kf_converged, "kf_converged/O");

  TChain* MC = NULL, *event = NULL, *jets = NULL, *MET = NULL, *muons = NULL, *electrons = NULL, *mtt = NULL;
  loadChain(inputFiles, MC, event, jets, MET, muons, electrons, mtt);

  const int nBins = 12;
  const double bins[] = {340, 400, 450, 500, 550, 600, 650, 750, 850, 950, 1050, 1200, 1400};

  GaussianProfile* h_mtt_reco_after_kf_vs_mtt = new GaussianProfile("mtt_reco_after_kf_vs_mtt", nBins, bins);
  h_mtt_reco_after_kf_vs_mtt->setPrefix("mtt");

  GaussianProfile* h_mtt_reco_vs_mtt = new GaussianProfile("mtt_reco_vs_mtt", nBins, bins);
  h_mtt_reco_vs_mtt->setPrefix("mtt");

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

  TH1F* h_min_kf_chisquare = new TH1F("min_kf_chisquare", "Minimum KinFit Chi square", 200, 0, 100);
  TH1F* h_kf_chisquare = new TH1F("kf_chisquare", "KinFit Chi squares", 400, 0, 200);

  MC->SetBranchStatus("*", 0);
  event->SetBranchStatus("*", 0);
  jets->SetBranchStatus("*", 0);
  MET->SetBranchStatus("*", 0);
  muons->SetBranchStatus("*", 0);
  electrons->SetBranchStatus("*", 0);
  mtt->SetBranchStatus("*", 0);

  // Jets
  SetBranchAddress(jets, "n_jets", &n_jets);
  SetBranchAddress(jets, "jet_4vector", &m_jet_lorentzvector);
  SetBranchAddress(jets, "jet_mcParticleIndex", &m_jet_MCIndex);

  // MTT
  SetBranchAddress(mtt, "MC_leptonicBIndex", &genLeptonicBIndex);
  SetBranchAddress(mtt, "MC_hadronicBIndex", &genHadronicBIndex);
  SetBranchAddress(mtt, "MC_hadronicFirstJetIndex", &genFirstJetIndex);
  SetBranchAddress(mtt, "MC_hadronicSecondJetIndex", &genSecondJetIndex);

  int selectedLeptonicBIndex, selectedHadronicBIndex, selectedFirstJetIndex, selectedSecondJetIndex;
  SetBranchAddress(mtt, "selectedLeptonicBIndex", &selectedLeptonicBIndex);
  SetBranchAddress(mtt, "selectedHadronicBIndex", &selectedHadronicBIndex);
  SetBranchAddress(mtt, "selectedHadronicFirstJetIndex", &selectedFirstJetIndex);
  SetBranchAddress(mtt, "selectedHadronicSecondJetIndex", &selectedSecondJetIndex);

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

  uint64_t selectedEntries = 0;

  uint64_t entries = MC->GetEntries();
  //uint64_t entries = 5000;
  if (ENTRIES > 0 && ENTRIES <= entries)
    entries = ENTRIES;

  // For efficiencies
  uint64_t numberOfSelectedEntries = 0;
  uint64_t numberOfAssociableEntries = 0;
  uint64_t numberOfMatchedEntries = 0;
  uint64_t numberOfWellPlacedEntries = 0;
  uint64_t numberOfConvergedEntries = 0;

  for (uint64_t entry = 0; entry < entries; entry++) {
    MC->GetEntry(entry);
    event->GetEntry(entry);
    jets->GetEntry(entry);
    MET->GetEntry(entry);
    muons->GetEntry(entry);
    electrons->GetEntry(entry);
    mtt->GetEntry(entry);

    if (((entry) % 500) == 0) {
      std::cout << "Processing entry " << entry + 1 << " out of " << entries << std::endl;
    }
    /*if (((entry) % 1) == 0) {
      std::cout << "Processing entry " << entry + 1 << " out of " << entries << std::endl;
    }*/

    if (numComb < 1)
      continue;
    
    // MET
    TLorentzVector* neutrinoP4 = static_cast<TLorentzVector*>((*m_met_lorentzvector)[0]);
    if (neutrinoP4->Pt() < 20)
      continue;

    // Lepton
    TLorentzVector* leptonP4 = NULL;
    int leptonCharge = 0;
    float ptLeptonCut = 0;
    if (SEMI_MU) {
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
    firstJetP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[selectedFirstJetIndex]);
    secondJetP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[selectedSecondJetIndex]);
    leptonicBP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[selectedLeptonicBIndex]);
    hadronicBP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[selectedHadronicBIndex]);

    int recoLeptonicBIndex = -1, recoHadronicBIndex = -1, recoFirstJetIndex = -1, recoSecondJetIndex = -1;

    // Is this solution the good one?
    if (!isMatched(genLeptonicBIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets) ||
        !isMatched(genHadronicBIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets) ||
        !isMatched(genFirstJetIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets)  ||
        !isMatched(genSecondJetIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets)) {
      goodCombinaison = false;
      associable = false;
    } else {
      // We know taht every partons have a matched jets.
      for (int i = 0; i < n_jets; i++) {
        int index = m_jet_MCIndex[i];
        if (index == genLeptonicBIndex)
          recoLeptonicBIndex = i;
        else if (index == genHadronicBIndex)
          recoHadronicBIndex = i;
        else if (index == genFirstJetIndex)
          recoFirstJetIndex = i;
        else if (index == genSecondJetIndex)
          recoSecondJetIndex = i;
      }

      associable = true;
      goodCombinaison = (recoLeptonicBIndex != -1) && (recoHadronicBIndex != -1) && (recoFirstJetIndex != -1) && (recoSecondJetIndex != -1);

      goodCombinaison =  (recoLeptonicBIndex == selectedLeptonicBIndex || recoLeptonicBIndex == selectedHadronicBIndex ||
                          recoLeptonicBIndex == selectedFirstJetIndex  || recoLeptonicBIndex == selectedSecondJetIndex);
      
      goodCombinaison &= (recoHadronicBIndex == selectedLeptonicBIndex || recoHadronicBIndex == selectedHadronicBIndex ||
                          recoHadronicBIndex == selectedFirstJetIndex  || recoHadronicBIndex == selectedSecondJetIndex);

      goodCombinaison &= (recoFirstJetIndex == selectedLeptonicBIndex || recoFirstJetIndex == selectedHadronicBIndex ||
                          recoFirstJetIndex == selectedFirstJetIndex  || recoFirstJetIndex == selectedSecondJetIndex);

      goodCombinaison &= (recoSecondJetIndex == selectedLeptonicBIndex || recoSecondJetIndex == selectedHadronicBIndex ||
                          recoSecondJetIndex == selectedFirstJetIndex  || recoSecondJetIndex == selectedSecondJetIndex);
    }


    if (! kinFitter.PzNeutrino(*leptonP4, *neutrinoP4, *leptonicBP4))
      continue;

    double mtt = (*firstJetP4 + *secondJetP4 + *hadronicBP4 + *leptonicBP4 + *leptonP4 + *neutrinoP4).M();
    h_mtt_reco_vs_mtt->fill(MC_mtt, mtt);
    mtt_reco = mtt;

    std::vector<TLorentzVector*> jets = {firstJetP4, secondJetP4, hadronicBP4, leptonicBP4};
    double minChiSquare = 1e15;
    kf_converged = false;

    reco::Candidate::PolarLorentzVector firstJetP4_after_kf, secondJetP4_after_kf, leptonicBP4_after_kf, hadronicBP4_after_kf, leptonP4_after_kf, neutrinoP4_after_kf;
    int kf_firstJetIndex = -1, kf_secondJetIndex = -1, kf_leptonicBIndex = -1, kf_hadronicBIndex = -1;

    for (int j1 = 0; j1 < 4; j1++) {
      for (int j2 = 0; j2 < 4; j2++) {
        if (j1 == j2)
          continue;

        for (int j3 = 0; j3 < 4; j3++) {
          if (j3 == j1 || j3 == j2)
            continue;

          for (int j4 = j3 + 1; j4 < 4; j4++) {
            if (j4 == j1 || j4 == j2 || j4 == j3)
              continue;

            if (fitter.fit(*jets[j3], *jets[j4], *jets[j2], *jets[j1], *leptonP4, *neutrinoP4, leptonCharge, (! SEMI_MU) ? CovarianceMatrix::kElectron : CovarianceMatrix::kMuon) != 0)
              continue;

            double chi2 = fitter.fitS();
            h_kf_chisquare->Fill(chi2);

            if (chi2 < minChiSquare) {
              minChiSquare = chi2;
              kf_converged = true;

              firstJetP4_after_kf = fitter.fittedHadP().p4();
              secondJetP4_after_kf = fitter.fittedHadQ().p4();
              hadronicBP4_after_kf = fitter.fittedHadB().p4();
              leptonicBP4_after_kf = fitter.fittedLepB().p4();
              leptonP4_after_kf = fitter.fittedLepton().p4();
              neutrinoP4_after_kf = fitter.fittedNeutrino().p4();

              kf_proba = fitter.fitProb();
              kf_chisquare = chi2;

              kf_leptonicBIndex = j1;
              kf_hadronicBIndex = j2;
              kf_firstJetIndex = j3;
              kf_secondJetIndex = j4;
            }
          }
        }
      }
    }

    numberOfSelectedEntries++;
    if (associable) {
      numberOfAssociableEntries++;

      if (goodCombinaison)
        numberOfMatchedEntries++;
    }

    if (! kf_converged) {
      tree->Fill();
      continue;
    }

    numberOfConvergedEntries++;

    // Look if the KinFit choose the combinaison where jets are well placed
    // HINT: B-tagging might help here
    // We have: std::vector<TLorentzVector*> jets = {firstJetP4, secondJetP4, hadronicBP4, leptonicBP4};
    // So well placed solution is: kf_leptonicBIndex == 3, kf_hadronicBIndex == 2, kf_firstJetIndex == 1/0, kf_secondJetIndex == 0/1 
    wellPlacedAfterKF = goodCombinaison && ((kf_leptonicBIndex == 3) && (kf_hadronicBIndex == 2) && (kf_firstJetIndex == 1 || kf_firstJetIndex == 0) && (kf_secondJetIndex == 0 || kf_secondJetIndex == 1));

    if (associable && goodCombinaison && wellPlacedAfterKF)
      numberOfWellPlacedEntries++;

    h_min_kf_chisquare->Fill(minChiSquare);

    selectedEntries++;

    double mtt_after_kf = (firstJetP4_after_kf + secondJetP4_after_kf + hadronicBP4_after_kf + leptonicBP4_after_kf + leptonicBP4_after_kf + neutrinoP4_after_kf).M();
    h_mtt_reco_after_kf_vs_mtt->fill(MC_mtt, mtt_after_kf);
    mtt_reco_after_kf = mtt_after_kf;


    // Delta
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

    tree->Fill();
  }

  std::cout << "Selected entries for mtt computation: " << selectedEntries << std::endl;

  delete MC; delete event; delete MET; delete jets; delete muons; delete electrons;

  if (outputFile.length() > 0) {
    TFile* f = TFile::Open(outputFile.c_str(), "recreate");
    h_mtt_reco_vs_mtt->write(f);
    h_mtt_reco_after_kf_vs_mtt->write(f);

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

    h_kf_chisquare->Write();
    h_min_kf_chisquare->Write();

    tree->Write();

    f->Close();
    delete f;
  }

  // Some efficiencies
  std::cout << "Associable efficiency: " << std::endl;
  printEff((double) numberOfAssociableEntries / numberOfSelectedEntries, numberOfSelectedEntries);

  std::cout << "Matched efficiency: " << std::endl;
  printEff((double) numberOfMatchedEntries / numberOfAssociableEntries, numberOfAssociableEntries);

  std::cout << "Well placed efficiency: " << std::endl;
  printEff((double) numberOfWellPlacedEntries / numberOfMatchedEntries, numberOfMatchedEntries);

  std::cout << "KinFit efficiency: " << std::endl;
  printEff((double) numberOfConvergedEntries / numberOfSelectedEntries, numberOfSelectedEntries);

  delete h_mtt_reco_vs_mtt;
  delete h_mtt_reco_after_kf_vs_mtt;

  delete h_delta_leptonicB_Et_vs_mtt;
  delete h_delta_leptonicB_Eta_vs_mtt;
  delete h_delta_leptonicB_Phi_vs_mtt;

  delete h_delta_hadronicB_Et_vs_mtt;
  delete h_delta_hadronicB_Eta_vs_mtt;
  delete h_delta_hadronicB_Phi_vs_mtt;

  delete h_delta_firstJet_Et_vs_mtt;
  delete h_delta_firstJet_Eta_vs_mtt;
  delete h_delta_firstJet_Phi_vs_mtt;

  delete h_delta_secondJet_Et_vs_mtt;
  delete h_delta_secondJet_Eta_vs_mtt;
  delete h_delta_secondJet_Phi_vs_mtt;

  delete h_delta_lepton_Et_vs_mtt;
  delete h_delta_lepton_Eta_vs_mtt;
  delete h_delta_lepton_Phi_vs_mtt;

  delete h_delta_neutrino_Et_vs_mtt;
  delete h_delta_neutrino_Eta_vs_mtt;
  delete h_delta_neutrino_Phi_vs_mtt;

  delete h_min_kf_chisquare;
  delete h_kf_chisquare;

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

    TCLAP::SwitchArg semimuArg("", "semi-mu", "Do analysis for semi-mu events", cmd);

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
    SEMI_MU  = semimuArg.getValue();
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
