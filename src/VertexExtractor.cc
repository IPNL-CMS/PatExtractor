#include "../interface/VertexExtractor.h"


VertexExtractor::VertexExtractor(bool doTree, edm::InputTag tag)
{
  m_tag = tag;

  // Set everything to 0
  VertexExtractor::reset();

  // Tree definition

  if (doTree)
  {
    m_tree_vtx      = new TTree("Vertices","RECO PV info") ;
    m_tree_vtx->Branch("n_vertices",      &m_n_vertices,   "n_vertices/I");  
    m_tree_vtx->Branch("vertex_vx",       &m_vtx_vx,       "vertex_vx[n_vertices]/F");  
    m_tree_vtx->Branch("vertex_vy",       &m_vtx_vy,       "vertex_vy[n_vertices]/F");  
    m_tree_vtx->Branch("vertex_vz",       &m_vtx_vz,       "vertex_vz[n_vertices]/F"); 
    m_tree_vtx->Branch("vertex_isFake",   &m_vtx_isFake,   "vertex_isFake[n_vertices]/B"); 
    m_tree_vtx->Branch("vertex_ndof",     &m_vtx_ndof,     "vertex_ndof[n_vertices]/F"); 
    m_tree_vtx->Branch("vertex_normChi2", &m_vtx_normChi2, "vertex_normChi2[n_vertices]/F");
    m_tree_vtx->Branch("vertex_ntracks",  &m_vtx_ntracks,  "vertex_ntracks[n_vertices]/I");
  }
}


VertexExtractor::VertexExtractor(TFile *a_file)
{
  std::cout << "VertexExtractor objet is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_vtx = dynamic_cast<TTree*>(a_file->Get("Vertices"));

  if (!m_tree_vtx)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_tree_vtx->SetBranchAddress("n_vertices",       &m_n_vertices);
  m_tree_vtx->SetBranchAddress("vertex_vx",        &m_vtx_vx);
  m_tree_vtx->SetBranchAddress("vertex_vy",        &m_vtx_vy);
  m_tree_vtx->SetBranchAddress("vertex_vz",        &m_vtx_vz);
  m_tree_vtx->SetBranchAddress("vertex_isFake",    &m_vtx_isFake);
  m_tree_vtx->SetBranchAddress("vertex_ndof",      &m_vtx_ndof);
  m_tree_vtx->SetBranchAddress("vertex_normChi2",  &m_vtx_normChi2);
  m_tree_vtx->SetBranchAddress("vertex_ntracks",   &m_vtx_ntracks);
}

VertexExtractor::~VertexExtractor()
{}



//
// Method filling the main particle tree
//

void VertexExtractor::writeInfo(const edm::Event *event) 
{
  edm::Handle<reco::VertexCollection> vertexHandle;
  event->getByLabel(m_tag, vertexHandle);
  const reco::VertexCollection vertexCollection = *(vertexHandle.product());

  VertexExtractor::reset();
  VertexExtractor::fillSize(static_cast<int>(vertexCollection.size()));

  if (VertexExtractor::getSize())
  {
    for(int i=0; i<VertexExtractor::getSize(); ++i) 
      VertexExtractor::writeInfo(&vertexCollection.at(i),i); 
  }

  VertexExtractor::fillTree();
}

void VertexExtractor::writeInfo(const reco::Vertex *part, int index) 
{
  if (index>=m_vertices_MAX) return;

  m_vtx_vx[index]      = part->position().x();
  m_vtx_vy[index]      = part->position().y();
  m_vtx_vz[index]      = part->position().z();
  m_vtx_isFake[index]  = part->isFake();
  m_vtx_ndof[index]    = part->ndof();
  m_vtx_normChi2[index]= part->normalizedChi2();
  m_vtx_ntracks[index] = static_cast<int>(part->tracksSize());
}


//
// Method getting the info from an input file
//

void VertexExtractor::getInfo(int ievt) 
{
  m_tree_vtx->GetEntry(ievt); 
}

// Method initializing everything (to do for each event)

void VertexExtractor::reset()
{
  m_n_vertices = 0;

  for (int i=0;i<m_vertices_MAX;++i) 
  {
    m_vtx_vx[i]        = 0.;
    m_vtx_vy[i]        = 0.;
    m_vtx_vz[i]        = 0.;
    m_vtx_isFake[i]    = 0;
    m_vtx_ndof[i]      = 0.;
    m_vtx_normChi2[i]  = 0.;
    m_vtx_ntracks[i]   = 0;
  }
}


void VertexExtractor::fillTree()
{
  m_tree_vtx->Fill(); 
}
 
void VertexExtractor::fillSize(int size)
{
  m_n_vertices=size;
}

int  VertexExtractor::getSize()
{
  return m_n_vertices;
}

float VertexExtractor::dist_to_vtx(int vtxidx, float x, float y, float z)
{
  float dx = vx(vtxidx)-x;
  float dy = vy(vtxidx)-y;
  float dz = vz(vtxidx)-z;

  return sqrt(dx*dx+dy*dy+dz*dz);

}
