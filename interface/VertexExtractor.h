#ifndef VERTEXEXTRACTOR_H
#define VERTEXEXTRACTOR_H

/**
 * VertexExtractor
 * \brief: Base class for extracting vertex info
 */


//Include RECO info
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class VertexExtractor
{

 public:

  VertexExtractor(bool doTree,edm::InputTag tag);
  VertexExtractor(TFile *a_file);
  ~VertexExtractor();


  void writeInfo(const reco::Vertex *part, int index); 
  void writeInfo(const edm::Event *event); 
  void getInfo(int ievt); 

  void reset();
  void fillTree(); 
  void fillSize(int size);
  float dist_to_vtx(int vtxidx, float x, float y, float z);

  int   getSize();

  float vx(int vtxidx)           {return m_vtx_vx[vtxidx];}
  float vy(int vtxidx)           {return m_vtx_vy[vtxidx];}
  float vz(int vtxidx)           {return m_vtx_vz[vtxidx];}
  float getVtxNdof(int vtxidx)   {return m_vtx_ndof[vtxidx];}
  bool  getVtxIsFake(int vtxidx) {return m_vtx_isFake[vtxidx];}
  int   getNtracks(int vtxidx)   {return m_vtx_ntracks[vtxidx];}
  float getNormChi2(int vtxidx)  {return m_vtx_normChi2[vtxidx];}

  bool isOK() {return m_OK;}

 private:
  
  TTree* m_tree_vtx;

  static const int 	m_vertices_MAX   = 50;

  bool m_OK;
  edm::InputTag m_tag;

  int   m_n_vertices;
  float	m_vtx_vx[m_vertices_MAX];
  float	m_vtx_vy[m_vertices_MAX];
  float	m_vtx_vz[m_vertices_MAX];
  bool  m_vtx_isFake[m_vertices_MAX];
  float m_vtx_ndof[m_vertices_MAX];
  float m_vtx_normChi2[m_vertices_MAX];
  int   m_vtx_ntracks[m_vertices_MAX];
};

#endif 
