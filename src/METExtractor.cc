#include "../interface/METExtractor.h"


METExtractor::METExtractor(bool doTree,edm::InputTag tag)
{
  // Set everything to 0

  m_tag = tag;
  m_met_lorentzvector = new TClonesArray("TLorentzVector");
  METExtractor::reset();

  // Tree definition

  if (doTree)
  {
    m_tree_met      = new TTree("MET_PF","PAT PF MET info");  
    m_tree_met->Branch("n_mets",  &m_n_mets,   "n_mets/I");  
    m_tree_met->Branch("met_4vector","TClonesArray",&m_met_lorentzvector, 1000, 0);
  }
}

METExtractor::METExtractor(TFile *a_file)
{
  std::cout << "METExtractor objet is retrieved" << std::endl;


  // Tree definition
  m_OK = false;

  m_tree_met = dynamic_cast<TTree*>(a_file->Get("MET_PF"));


  if (!m_tree_met)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_met_lorentzvector = new TClonesArray("TLorentzVector");

  if (m_tree_met->FindBranch("n_mets"))
    m_tree_met->SetBranchAddress("n_mets",  &m_n_mets);

  if (m_tree_met->FindBranch("met_4vector"))
    m_tree_met->SetBranchAddress("met_4vector",&m_met_lorentzvector);

}

METExtractor::~METExtractor()
{}



//
// Method filling the main particle tree
//

void METExtractor::writeInfo(const edm::Event *event) 
{
  edm::Handle< edm::View<pat::MET> >  METHandle;
  event->getByLabel(m_tag, METHandle);
  edm::View<pat::MET> p_METs = *METHandle;

  METExtractor::reset();
  METExtractor::fillSize(static_cast<int>(p_METs.size()));

  if (METExtractor::getSize())
  {
    for(int i=0; i<METExtractor::getSize(); ++i) 
      METExtractor::writeInfo(&p_METs.at(i),i); 
  }

  METExtractor::fillTree();
}


void METExtractor::writeInfo(const pat::MET *part, int index) 
{
  if (index>=m_mets_MAX) return;

  new((*m_met_lorentzvector)[index]) TLorentzVector(part->px(),part->py(),part->pz(),part->energy());
}


//
// Method getting the info from an input file
//

void METExtractor::getInfo(int ievt) 
{
  m_tree_met->GetEntry(ievt); 
}

// Method initializing everything (to do for each event)

void METExtractor::reset()
{
  m_n_mets = 0;

  m_met_lorentzvector->Clear();
}


void METExtractor::fillTree()
{
  m_tree_met->Fill(); 
}
 
void METExtractor::fillSize(int size)
{
  m_n_mets=size;
}

int  METExtractor::getSize()
{
  return m_n_mets;
}

void METExtractor::setMETLorentzVector(int idx, float E, float Px, float Py, float Pz)
{
  new((*m_met_lorentzvector)[idx]) TLorentzVector(Px,Py,Pz,E);
}
