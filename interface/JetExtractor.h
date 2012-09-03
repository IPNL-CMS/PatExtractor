#ifndef JETEXTRACTOR_H
#define JETEXTRACTOR_H

/**
 * JetExtractor
 * \brief: Base class for extracting jet info
 */


//Include PAT info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <CommonTools/Utils/interface/PtComparator.h>

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/View.h"

#include "../interface/MCExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class JetExtractor
{

 public:

  JetExtractor(bool doTree, edm::InputTag tag, bool correctJets, const std::string& jetCorrectorLabel);
  JetExtractor(TFile *a_file);
  ~JetExtractor();


  void writeInfo(const edm::Event *event, const edm::EventSetup& iSetup, MCExtractor* m_MC, bool doMC); 

  void writeInfo(const pat::Jet *part, int index); 
  void getInfo(int ievt); 

  void reset();
  void fillTree(); 
  void fillSize(int size);
  int  getSize();

  int getMatch(const pat::Jet *part, MCExtractor* m_MC);

  TLorentzVector *getJetLorentzVector(int jetidx) {return (TLorentzVector*)m_jet_lorentzvector->At(jetidx);}

  float getJetBTagProb_CSV(int jetidx) const { return m_jet_btag_CSV[jetidx]; }
  float getJetBTagProb_TCHP(int jetidx) const { return m_jet_btag_TCHP[jetidx]; }

  //float getJetBTagProb_SSVHP(int jetidx) {return m_jet_btag_SSVHP[jetidx];}
  //float getJetBTagProb_TCHE(int jetidx) {return m_jet_btag_TCHE[jetidx];}
  
  int getJetMCIndex(int jetidx){return m_jet_MCIndex[jetidx];}
  
  bool isOK() {return m_OK;}

  /**
   * This method is need for inline JEC. We need to be able to change jet pt on-the-fly
   */
  void setJetLorentzVector(int index, TLorentzVector* vector)
  {
    if (index >= m_jet_lorentzvector->GetEntries())
      return;

    new((*m_jet_lorentzvector)[index]) TLorentzVector(*vector);
  }

  void setJetLorentzVector(int index, double energy, double px, double py, double pz)
  {
    TLorentzVector* v = new TLorentzVector(px, py, pz, energy);
    setJetLorentzVector(index, v);
    delete v;
  }

  // Jet ID
  bool isPFJetLoose(const pat::Jet& jet);

  pat::JetCollection correctJets(pat::JetCollection& jets, const edm::Event& iEvent, const edm::EventSetup& iSetup);

 private:
  
  TTree* m_tree_jet;

  bool mCorrectJets;
  std::string mJetCorrectorLabel;
  GreaterByPt<pat::Jet> mSorter;

  static const int 	m_jets_MAX       = 200;

  edm::InputTag m_tag;
  float m_deltaR_cut;

  bool m_OK;
  int   m_n_jets;

  TClonesArray* m_jet_lorentzvector;
  float	m_jet_vx[m_jets_MAX];
  float	m_jet_vy[m_jets_MAX];
  float	m_jet_vz[m_jets_MAX];
  int	m_jet_chmult[m_jets_MAX];
  float	m_jet_chmuEfrac[m_jets_MAX];
  float	m_jet_chemEfrac[m_jets_MAX];
  float	m_jet_chhadEfrac[m_jets_MAX];
  float	m_jet_nemEfrac[m_jets_MAX];
  float	m_jet_nhadEfrac[m_jets_MAX];

  //float	m_jet_btag_BjetProb[m_jets_MAX];
  //float	m_jet_btag_SSVHE[m_jets_MAX];
  //float	m_jet_btag_SSVHP[m_jets_MAX];
  //float	m_jet_btag_TCHE[m_jets_MAX];
  float	m_jet_btag_jetProb[m_jets_MAX];
  float	m_jet_btag_TCHP[m_jets_MAX];
  float	m_jet_btag_CSV[m_jets_MAX];

  int   m_jet_MCIndex[m_jets_MAX];
};

#endif 
