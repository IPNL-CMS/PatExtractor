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

#include "Extractors/PatExtractor/interface/BaseExtractor.h"

//Include std C++
#include <iostream>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

class VertexExtractor: public BaseExtractor<reco::Vertex>
{

  public:

    VertexExtractor(const std::string& name, const edm::InputTag& tag, bool doTree);
    VertexExtractor(const std::string& name, TFile *a_file);
    virtual ~VertexExtractor();

    void writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const reco::Vertex& part, int index); 
    void getInfo(int ievt); 

    virtual const reco::Candidate* getGenParticle(const reco::Vertex& vtx) {
      return NULL;
    }

    virtual void setGenParticleIndex(int genParticleIndex, int index) {
    }

    virtual float getMCMatchDeltaR() {
      return 0.;
    }

    virtual float getMCMatchDPtRel() {
      return 0.;
    }

    virtual std::vector<int> getPdgIds() {
      return {};
    }

    virtual TLorentzVector getP4(const reco::Vertex& object) {
      TLorentzVector p4;
      return p4;
    }

    void reset();
    void fillTree(); 
    float dist_to_vtx(int vtxidx, float x, float y, float z);

    float vx(int vtxidx)           {return m_vtx_vx[vtxidx];}
    float vy(int vtxidx)           {return m_vtx_vy[vtxidx];}
    float vz(int vtxidx)           {return m_vtx_vz[vtxidx];}
    float getVtxNdof(int vtxidx)   {return m_vtx_ndof[vtxidx];}
    bool  getVtxIsFake(int vtxidx) {return m_vtx_isFake[vtxidx];}
    int   getNtracks(int vtxidx)   {return m_vtx_ntracks[vtxidx];}
    float getNormChi2(int vtxidx)  {return m_vtx_normChi2[vtxidx];}

  private:

    TTree* m_tree_vtx;

    static const int 	m_vertices_MAX   = 50;

    float	m_vtx_vx[m_vertices_MAX];
    float	m_vtx_vy[m_vertices_MAX];
    float	m_vtx_vz[m_vertices_MAX];
    bool  m_vtx_isFake[m_vertices_MAX];
    float m_vtx_ndof[m_vertices_MAX];
    float m_vtx_normChi2[m_vertices_MAX];
    int   m_vtx_ntracks[m_vertices_MAX];
};

#endif 
