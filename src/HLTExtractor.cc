#include "../interface/HLTExtractor.h"


HLTExtractor::HLTExtractor(bool doTree)
{
  // Set everything to 0

  m_OK         = false;
  m_HLT_vector = new std::vector< std::string >;
  HLTExtractor::reset();

  // Tree definition

  if (doTree)
  {
    m_OK = true;
    m_tree_HLT       = new TTree("HLT","HLT info");  
    m_tree_HLT->Branch("n_paths",  &m_n_HLTs,"n_paths/I");       
    m_tree_HLT->Branch("HLT_vector","vector<string>",&m_HLT_vector);
  }

}

HLTExtractor::HLTExtractor(TFile *a_file)
{
  std::cout << "HLTExtractor objet is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_HLT = dynamic_cast<TTree*>(a_file->Get("HLT"));

  if (!m_tree_HLT)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  // Branches definition

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

void HLTExtractor::writeInfo(const edm::Event *event) 
{
  HLTExtractor::reset();

  edm::Handle<edm::TriggerResults> triggerResults ;
  edm::InputTag tag("TriggerResults", "", "HLT");
  event->getByLabel(tag,triggerResults);

  if (triggerResults.isValid())
  {
    const edm::TriggerNames & triggerNames = event->triggerNames(*triggerResults);
    
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


  HLTExtractor::fillTree();
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
 
void HLTExtractor::fillSize(int size)
{
  m_n_HLTs=size;
}

int  HLTExtractor::getSize()
{
  return m_n_HLTs;
}
 
 
void HLTExtractor::print()
{
  if (HLTExtractor::isOK())
  {
    std::cout << "------------------------------------" << std::endl;
    std::cout << "This event fired " << HLTExtractor::getSize() << " HLT path(s):" << std::endl;
    
    for(int i=0 ; i<HLTExtractor::getSize(); ++i) 
    {
      std::cout << i+1 << " : " << HLTExtractor::paths(i) << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;
  }
  else
  {
    std::cout << "There is no HLT tree in this rootuple" << std::endl;
  }
}
