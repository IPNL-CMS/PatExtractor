#include "../interface/EventExtractor.h"

EventExtractor::EventExtractor(const std::string& name, const edm::ParameterSet& parameters):
  SuperBaseExtractor(name, parameters),
  m_puSummaryTag(parameters.getParameter<edm::InputTag>("pileup_summary")),
  m_generatorTag(parameters.getParameter<edm::InputTag>("generator"))
{
  // Tree definition
  m_OK = true;

  m_tree_event    = new TTree(name.c_str(), "Event info");  

  // Branches definition

  m_tree_event->Branch("evtID",  &m_evtID,"evtID/i");
  m_tree_event->Branch("lumi",   &m_lumi,"lumi/i");
  m_tree_event->Branch("run",    &m_run,"run/i");
  m_tree_event->Branch("BCID",   &m_BCID,"BCID/I");
  m_tree_event->Branch("time",   &m_time,"time/i");
  m_tree_event->Branch("nPU",    &m_nPU,"nPileUp/I");
  m_tree_event->Branch("nTrueInteractions", &m_nTrueInteractions, "nTrueInteractions/F");

  m_tree_event->Branch("generator_weight", &m_generator_weight, "generator_weight/F");

  // Set everything to 0

  EventExtractor::reset();
}

EventExtractor::EventExtractor(const std::string& name, const edm::ParameterSet& settings, TFile *file):
  SuperBaseExtractor(name, settings, file)
{
  // Tree definition

  m_tree_event = dynamic_cast<TTree*>(file->Get(name.c_str()));

  m_OK = false;

  if (!m_tree_event)
    std::cout << "Event tree not defined, this is bad" << std::endl;

  // Branches definition

  m_n_events = m_tree_event->GetEntries();

  std::cout << "This file contains " << m_n_events << " events..." << std::endl;

  if (m_tree_event->FindBranch("evtID"))   m_tree_event->SetBranchAddress("evtID",  &m_evtID);
  if (m_tree_event->FindBranch("lumi"))    m_tree_event->SetBranchAddress("lumi",   &m_lumi);
  if (m_tree_event->FindBranch("run"))     m_tree_event->SetBranchAddress("run",    &m_run);
  if (m_tree_event->FindBranch("BCID"))    m_tree_event->SetBranchAddress("BCID",   &m_BCID);
  if (m_tree_event->FindBranch("time"))    m_tree_event->SetBranchAddress("time",   &m_time);
  if (m_tree_event->FindBranch("nPileUp")) m_tree_event->SetBranchAddress("nPileUp",&m_nPU);
  if (m_tree_event->FindBranch("nTrueInteractions"))
    m_tree_event->SetBranchAddress("nTrueInteractions", &m_nTrueInteractions);
  if (m_tree_event->FindBranch("generator_weight"))
    m_tree_event->SetBranchAddress("generator_weight", &m_generator_weight);

  m_OK = true;
}

EventExtractor::~EventExtractor()
{}

void EventExtractor::doConsumes(edm::ConsumesCollector&& collector) {
  SuperBaseExtractor::doConsumes(std::forward<edm::ConsumesCollector>(collector));

  if (m_isMC) {
    m_puSummaryToken = collector.consumes<std::vector<PileupSummaryInfo>>(m_puSummaryTag);
    m_generatorToken = collector.consumes<GenEventInfoProduct>(m_generatorTag);
  }
}

//
// Method filling the main particle tree
//

void EventExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor)
{
  m_evtID             = (event.eventAuxiliary()).id().event();
  m_BCID              = (event.eventAuxiliary()).bunchCrossing();
  m_time              = (event.eventAuxiliary()).time().unixTime();
  m_lumi              = (event.eventAuxiliary()).luminosityBlock();
  m_run               = (event.eventAuxiliary()).run();
  m_nPU               = -1;
  m_nTrueInteractions = -1;

  if (mcExtractor) 
  {
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    event.getByToken(m_puSummaryToken, PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();
      if(BX == 0) { 
        m_nPU = PVI->getPU_NumInteractions();
        m_nTrueInteractions = PVI->getTrueNumInteractions();
        continue;
      }
    }
  }

  if (m_isMC) {
    edm::Handle<GenEventInfoProduct> generatorInfo;
    event.getByToken(m_generatorToken, generatorInfo);

    if (generatorInfo.isValid()) {
      m_generator_weight = generatorInfo->weight();
    }
  }

  m_tree_event->Fill();
}


//
// Method getting the info from an input file
//

void EventExtractor::getInfo(int ievt) 
{
  m_tree_event->GetEntry(ievt); 
}

// Method initializing everything (to do for each event)

void EventExtractor::reset()
{
  m_evtID             =  0;
  m_BCID              =  0;
  m_time              =  0;
  m_lumi              =  0;
  m_run               =  0;
  m_nPU               =  0;
  m_nTrueInteractions = 0;
  m_generator_weight  = 1;
}

// Method print the event info

void EventExtractor::print()
{
  std::cout << "------------------------------------" << std::endl;
  std::cout << "Dealing with event : " << EventExtractor::evtID() << std::endl;
  std::cout << "From run           : " << EventExtractor::run() << std::endl;
  std::cout << "------------------------------------" << std::endl;
}

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, EventExtractor, "event_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, EventExtractor, "event_extractor");
