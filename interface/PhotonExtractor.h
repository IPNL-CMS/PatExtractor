#ifndef PHOTONEXTRACTOR_H
#define PHOTONEXTRACTOR_H

/**
 * PhotonExtractor
 * \brief: Base class for extracting photon info
 */


//Include PAT info
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "../interface/BaseExtractor.h"
#include "../interface/MCExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class PhotonExtractor: public BaseExtractor<pat::Photon>
{

 public:

  PhotonExtractor(const std::string& name, const edm::InputTag& tag, bool doTree);
  PhotonExtractor(const std::string& name, TFile *a_file);
  virtual ~PhotonExtractor();

  void writeInfo(const pat::Photon& part, int index); 
  void getInfo(int ievt);

  void reset();
  void fillTree(); 
  
  virtual void doMCMatch(const pat::Photon& object, MCExtractor* mcExtractor, int index);

  // Setters/Getters

 private:
  
  TTree* m_tree_photon;

  static const int 	m_photons_MAX    = 100;

  float m_deltaR_cut;

  TClonesArray* m_pho_lorentzvector;
  float	m_pho_vx[m_photons_MAX];
  float	m_pho_vy[m_photons_MAX];
  float	m_pho_vz[m_photons_MAX];
  int   m_pho_MCIndex[m_photons_MAX];
};

#endif 
