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

  TChain* MC = NULL, *event = NULL, *jets = NULL, *MET = NULL, *muons = NULL, *electrons = NULL, *mtt = NULL;
  loadChain(inputFiles, MC, event, jets, MET, muons, electrons, mtt);

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

  uint64_t possibleMatchableEntries = 0;

  uint64_t kf_totalCombinaison = 0;
  uint64_t kf_correctlyChoosedCombinaison = 0;
  std::vector<uint64_t> kf_goodCombinaisonPosition(24);

  //uint64_t entries = MC->GetEntries();
  uint64_t entries = 100000;

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

    std::vector<TLorentzVector*> jets = {firstJetP4, secondJetP4, hadronicBP4, leptonicBP4};
    std::vector<double> chiSquares;
    int goodCombinaisonIndex = -1;
    double minChiSquare = 1e15;

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

            if (fitter.fit(*jets[j3], *jets[j4], *jets[j2], *jets[j1], *leptonP4, *neutrinoP4, leptonCharge, (channel == 1) ? CovarianceMatrix::kElectron : CovarianceMatrix::kMuon) != 0)
              continue;

            double chi2 = fitter.fitS();
            chiSquares.push_back(chi2);

            //std::cout << "For [" << j1 << " " << j2 << " " << j3 << " " << j4 << "] chisquare: " << chi2 << std::endl;

            if ((j4 == 0 || j4 == 1) && (j3 == 0 || j3 == 1) && j2 == 2 && j1 == 3)
              goodCombinaisonIndex = chiSquares.size() - 1;

            if (chi2 < minChiSquare)
              minChiSquare = chi2;
          }
        }
      }
    }

    kf_totalCombinaison++;
    if (goodCombinaisonIndex >= 0 && (fabs(minChiSquare - chiSquares[goodCombinaisonIndex]) < 1e-6))
      kf_correctlyChoosedCombinaison++;

    if (goodCombinaisonIndex >= 0) {
      double goodChiSquare = chiSquares[goodCombinaisonIndex];
      std::sort(chiSquares.begin(), chiSquares.end());
      for (unsigned int i = 0; i < chiSquares.size(); i++) {
        if (fabs(goodChiSquare - chiSquares[i]) < 1e-6) {
          kf_goodCombinaisonPosition[i]++;
          break;
        }
      }
    }
  }

  std::cout << "Efficiency of choosing the right jet at the right place: " << (double) kf_correctlyChoosedCombinaison / kf_totalCombinaison * 100 << " % (" << kf_correctlyChoosedCombinaison << " / " << kf_totalCombinaison << ")" << std::endl;

  for (unsigned int i = 0; i < kf_goodCombinaisonPosition.size(); i++) {
    std::cout << "Number of time in #" << i << " position: " << (double) kf_goodCombinaisonPosition[i] / kf_totalCombinaison * 100 << " %" << std::endl;
  }

  delete MC; delete event; delete MET; delete jets; delete muons; delete electrons;

  if (outputFile.length() > 0) {
    TFile* f = TFile::Open(outputFile.c_str(), "recreate");
    f->Close();
    delete f;
  }

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
