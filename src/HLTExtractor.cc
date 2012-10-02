#include "../interface/HLTExtractor.h"


HLTExtractor::HLTExtractor(const std::string& name, bool doTree)
{
  // Set everything to 0

  m_OK         = false;
  m_HLT_vector = new std::vector< std::string >;
  reset();

  // Tree definition

  if (doTree)
  {
    m_OK = true;
    m_tree_HLT       = new TTree(name.c_str(), "HLT info");  
    m_tree_HLT->Branch("n_paths",  &m_n_HLTs,"n_paths/I");       
    m_tree_HLT->Branch("HLT_vector","vector<string>",&m_HLT_vector);
  }

}

HLTExtractor::HLTExtractor(const std::string& name, TFile *a_file)
{
  std::cout << "HLTExtractor objet is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_HLT = dynamic_cast<TTree*>(a_file->Get(name.c_str()));

  if (!m_tree_HLT)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  // Branches definition
  m_HLT_vector = new std::vector<std::string>();

  if (m_tree_HLT->FindBranch("n_paths")) 
    m_tree_HLT->SetBranchAddress("n_paths",  &m_n_HLTs);       

  if (m_tree_HLT->FindBranch("HLT_vector"))
    m_tree_HLT->SetBranchAddress("HLT_vector",&m_HLT_vector);
}

HLTExtractor::~HLTExtractor()
{}



//
// Method filling the main particle tree
//

void HLTExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor)
{
  reset();

  edm::Handle<edm::TriggerResults> triggerResults ;
  edm::InputTag tag("TriggerResults", "", "HLT");
  event.getByLabel(tag,triggerResults);

  if (triggerResults.isValid())
  {
    const edm::TriggerNames & triggerNames = event.triggerNames(*triggerResults);

    for(int i = 0 ; i < static_cast<int>(triggerResults->size()); i++) 
    {
      if (triggerResults->accept(i)!=0)
      {
        if (triggerNames.triggerName(i) == "HLTriggerFinalPath") continue; // This one is pretty useless...
        if ((triggerNames.triggerName(i).c_str())[0] == 'A') continue;     // Remove AlCa HLT paths

        m_HLT_vector->push_back(triggerNames.triggerName(i));
        m_n_HLTs++;
      }
    }
  }


  fillTree();
}


//
// Method getting the info from an input file
//

void HLTExtractor::getInfo(int ievt) 
{
  m_tree_HLT->GetEntry(ievt); 
}

// Method initializing everything (to do for each HLT)

void HLTExtractor::reset()
{
  m_n_HLTs = 0;
  m_HLT_vector->clear();
}


void HLTExtractor::fillTree()
{
  m_tree_HLT->Fill(); 
}

int HLTExtractor::getSize() const {
  return m_n_HLTs;
}
