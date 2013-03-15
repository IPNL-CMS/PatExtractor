#include "../interface/HLTExtractor.h"

#include "tinyxml2.h"

using namespace tinyxml2;

HLTExtractor::HLTExtractor(const std::string& name, bool doTree, const edm::ParameterSet& config):
  m_triggersXML(config.getUntrackedParameter<std::string>("triggersXML", ""))
{
  m_filterHLT = m_triggersXML.length() > 0;

  // Set everything to 0

  m_OK         = false;
  m_HLT_vector = new std::vector< std::string >;
  m_mustPass = new std::string();
  reset();

  // Tree definition

  if (doTree)
  {
    m_OK = true;
    m_tree_HLT       = new TTree(name.c_str(), "HLT info");  
    m_tree_HLT->Branch("n_paths",  &m_n_HLTs,"n_paths/I");       
    m_tree_HLT->Branch("HLT_vector","vector<string>",&m_HLT_vector);

    m_tree_HLT->Branch("HLT_filtered", &m_filterHLT, "HLT_filtered/O");
    m_tree_HLT->Branch("HLT_mustPass", &m_mustPass);
    m_tree_HLT->Branch("HLT_passed", &m_passed, "HLT_passed/O");
  }

  std::cout << std::endl;
  std::cout << "Triggers for this analysis:" << std::endl;
  if (m_filterHLT) {
    // Read XML document
    m_triggersService.reset(new Triggers(m_triggersXML));
    m_triggersService->print();
    std::cout << std::endl;
  } else {
    std::cout << "\tNo triggers" << std::endl;
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
  m_mustPass = NULL;

  if (m_tree_HLT->FindBranch("n_paths")) 
    m_tree_HLT->SetBranchAddress("n_paths",  &m_n_HLTs);       

  if (m_tree_HLT->FindBranch("HLT_vector"))
    m_tree_HLT->SetBranchAddress("HLT_vector",&m_HLT_vector);

  if (m_tree_HLT->FindBranch("HLT_filtered"))
    m_tree_HLT->SetBranchAddress("HLT_filtered", &m_filterHLT);

  if (m_tree_HLT->FindBranch("HLT_mustPass"))
    m_tree_HLT->SetBranchAddress("HLT_mustPass", &m_mustPass);

  if (m_tree_HLT->FindBranch("HLT_passed"))
    m_tree_HLT->SetBranchAddress("HLT_passed", &m_passed);
}

HLTExtractor::~HLTExtractor()
{}



//
// Method filling the main particle tree
//

void HLTExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor)
{
  reset();

  const boost::regex* triggerRegex = nullptr;
  if (m_filterHLT) {
    const PathVector& triggers = m_triggersService->getTriggers(event.run());
    *m_mustPass = triggers[0].str();
    triggerRegex = &triggers[0];
  }

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
        std::string triggerName = triggerNames.triggerName(i);
        if (triggerName == "HLTriggerFinalPath") continue; // This one is pretty useless...
        if (triggerName[0] == 'A') continue;     // Remove AlCa HLT paths

        m_HLT_vector->push_back(triggerName);
        m_n_HLTs++;

        if (triggerRegex) {
          if (boost::regex_match(triggerName, *triggerRegex)) {
            m_passed = true;
          }
        }
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
  m_passed = false;
  *m_mustPass = "";
}


void HLTExtractor::fillTree()
{
  m_tree_HLT->Fill(); 
}

int HLTExtractor::getSize() const {
  return m_n_HLTs;
}



bool Triggers::parse(const std::string& xmlContent) {
  XMLDocument doc;
  if (doc.Parse(xmlContent.c_str())) {
    doc.PrintError();
    return false;
  }

  const XMLElement* root = doc.FirstChildElement("triggers");
  if (! root)
    return false;

  const XMLElement* runs = root->FirstChildElement("runs");
  for (; runs; runs = runs->NextSiblingElement("runs")) {
    parseRunsElement(runs);
  }

  return true;
}

bool Triggers::parseRunsElement(const XMLElement* runs) {
  Range<unsigned> runRange(runs->UnsignedAttribute("from"), runs->UnsignedAttribute("to"));

  PathVector runPaths;

  const XMLElement* paths = runs->FirstChildElement("path");
  for (; paths; paths = paths->NextSiblingElement("path")) {
    const std::string name = paths->FirstChildElement("name")->GetText();
    runPaths.push_back(boost::regex(name, boost::regex_constants::icase));
  }

  mTriggers[runRange] = runPaths;
  return true;
}

const PathVector& Triggers::getTriggers(unsigned int run) {
  if (mCachedRange && mCachedRange->in(run)) {
    return *mCachedVector;
  }

  for (auto& trigger: mTriggers) {
    const Range<unsigned int>& runRange = trigger.first;

    if (runRange.in(run)) {

      mCachedRange = &runRange;
      mCachedVector = &trigger.second;

      return *mCachedVector;
    }
  }

  std::cout << "Error: run " << run << " not found for triggers selection" << std::endl;
  assert(false);
}

void Triggers::print() {
  for (auto& trigger: mTriggers) {
    const Range<unsigned int>& runRange = trigger.first;
    const auto& paths = trigger.second;

    std::cout << "Runs: " << runRange << std::endl;
    for (auto& path: paths) {
      std::cout << path << std::endl;
    }
  }
}
