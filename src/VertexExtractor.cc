#include "../interface/VertexExtractor.h"


VertexExtractor::VertexExtractor(const std::string& name, const edm::InputTag& tag, bool doTree)
  : BaseExtractor(name)
{
  m_tag = tag;

  // Set everything to 0
  reset();
  m_OK = false;

  // Tree definition

  if (doTree)
  {
    m_OK = true;
    m_tree_vtx      = new TTree(m_name.c_str(), "RECO PV info") ;
    m_tree_vtx->Branch("n_vertices",      &m_size,   "n_vertices/i");  
    m_tree_vtx->Branch("vertex_vx",       &m_vtx_vx,       "vertex_vx[n_vertices]/F");  
    m_tree_vtx->Branch("vertex_vy",       &m_vtx_vy,       "vertex_vy[n_vertices]/F");  
    m_tree_vtx->Branch("vertex_vz",       &m_vtx_vz,       "vertex_vz[n_vertices]/F"); 
    m_tree_vtx->Branch("vertex_isFake",   &m_vtx_isFake,   "vertex_isFake[n_vertices]/B"); 
    m_tree_vtx->Branch("vertex_ndof",     &m_vtx_ndof,     "vertex_ndof[n_vertices]/F"); 
    m_tree_vtx->Branch("vertex_normChi2", &m_vtx_normChi2, "vertex_normChi2[n_vertices]/F");
    m_tree_vtx->Branch("vertex_ntracks",  &m_vtx_ntracks,  "vertex_ntracks[n_vertices]/I");
  }
}


VertexExtractor::VertexExtractor(const std::string& name, TFile *a_file)
  : BaseExtractor(name)
{
  std::cout << "VertexExtractor objet is retrieved" << std::endl;
  m_file = a_file;

  // Tree definition
  m_OK = false;

  m_tree_vtx = dynamic_cast<TTree*>(a_file->Get(m_name.c_str()));

  if (!m_tree_vtx)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_tree_vtx->SetBranchAddress("n_vertices",       &m_size);
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

void VertexExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const reco::Vertex& part, int index) 
{
  if (index>=m_vertices_MAX) return;

  m_vtx_vx[index]      = part.position().x();
  m_vtx_vy[index]      = part.position().y();
  m_vtx_vz[index]      = part.position().z();
  m_vtx_isFake[index]  = part.isFake();
  m_vtx_ndof[index]    = part.ndof();
  m_vtx_normChi2[index]= part.normalizedChi2();
  m_vtx_ntracks[index] = static_cast<int>(part.tracksSize());
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
  m_size = 0;

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
 
float VertexExtractor::dist_to_vtx(int vtxidx, float x, float y, float z)
{
  float dx = vx(vtxidx)-x;
  float dy = vy(vtxidx)-y;
  float dz = vz(vtxidx)-z;

  return sqrt(dx*dx+dy*dy+dz*dz);

}
