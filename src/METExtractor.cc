#include "../interface/METExtractor.h"


METExtractor::METExtractor(const std::string& name, const edm::InputTag& tag, bool doTree)
  : BaseExtractor(name)
{
  // Set everything to 0

  m_tag = tag;
  m_met_lorentzvector = new TClonesArray("TLorentzVector");
  reset();

  // Tree definition
  m_OK = false;

  if (doTree)
  {
    m_OK = true;
    m_tree_met      = new TTree(name.c_str(), "PAT PF MET info");  
    m_tree_met->Branch("n_mets",  &m_size,   "n_mets/I");  
    m_tree_met->Branch("met_4vector","TClonesArray",&m_met_lorentzvector, 1000, 0);
  }
}

METExtractor::METExtractor(const std::string& name, TFile *a_file)
  : BaseExtractor(name)
{
  std::cout << "METExtractor objet is retrieved" << std::endl;
  m_file = a_file;

  // Tree definition
  m_OK = false;
  m_tree_met = dynamic_cast<TTree*>(a_file->Get(name.c_str()));

  if (!m_tree_met)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_met_lorentzvector = new TClonesArray("TLorentzVector");

  if (m_tree_met->FindBranch("n_mets"))
    m_tree_met->SetBranchAddress("n_mets",  &m_size);

  if (m_tree_met->FindBranch("met_4vector"))
    m_tree_met->SetBranchAddress("met_4vector",&m_met_lorentzvector);

}

METExtractor::~METExtractor()
{}

//
// Method filling the main particle tree
//

void METExtractor::writeInfo(const pat::MET& part, int index) 
{
  if (index>=m_mets_MAX) return;

  new((*m_met_lorentzvector)[index]) TLorentzVector(part.px(),part.py(),part.pz(),part.energy());
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
  m_size = 0;

  m_met_lorentzvector->Clear();
}


void METExtractor::fillTree()
{
  m_tree_met->Fill(); 
}
 
void METExtractor::setMETLorentzVector(int idx, float E, float Px, float Py, float Pz)
{
  new((*m_met_lorentzvector)[idx]) TLorentzVector(Px,Py,Pz,E);
}
