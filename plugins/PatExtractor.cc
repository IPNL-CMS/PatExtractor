#include "../interface/PatExtractor.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;
using namespace edm;

PatExtractor::PatExtractor(const edm::ParameterSet& config) :
  is_MC_         (config.getUntrackedParameter<bool>("isMC", true)),
  do_fill_       (config.getUntrackedParameter<bool>("fillTree", true)),
  nevts_         (config.getUntrackedParameter<int>("n_events", 10000)),
  outFilename_   (config.getParameter<std::string>("extractedRootFile")),
  inFilename_    (config.getParameter<std::string>("inputRootFile"))
{
  // Initializations

  nevent_tot = 0;

  // If do_fill is set to True, you extract the whole data, otherwise you start 
  // from a file already extracted (inFilename_)

  ScaleFactorService::createInstance(config);

  (do_fill_) 
    ? PatExtractor::initialize(config)
    : PatExtractor::retrieve(config);

  // Load plugins
  if (config.existsAs<edm::ParameterSet>("plugins")) {
    const edm::ParameterSet& plugins = config.getParameterSet("plugins");
    std::vector<std::string> pluginNames = plugins.getParameterNames();
    for (std::string& pluginName: pluginNames) {
      edm::ParameterSet pluginParameters;
      if (plugins.existsAs<edm::ParameterSet>(pluginName))
        pluginParameters = plugins.getParameterSet(pluginName);

      m_plugins.push_back(std::shared_ptr<patextractor::Plugin>(PatExtractorPluginFactory::get()->create(pluginName, pluginParameters)));
      m_plugins.back()->setIsMC(is_MC_);
    }
  }
}


// Job initialization
void PatExtractor::beginJob()
{
  for (auto& extractor: m_extractors) {
    extractor->beginJob(!do_fill_);
  }
}


// What to do at the start of a new run
void PatExtractor::beginRun(Run const& run, EventSetup const& setup)
{

  if (is_MC_)
    ScaleFactorService::getInstance().prepareBTaggingScaleFactors(setup);

  nevent = 0;

  // If we start from existing file we don't have to loop over events
  std::shared_ptr<EventExtractor> eventExtractor = std::static_pointer_cast<EventExtractor>(getExtractor("event"));
  if (!do_fill_ && eventExtractor->n_events()) 
  {    
    // If you start from an extracted file, the number of events you want to loop on
    // is defined as an option, not in CMSSW...

    nevent = min(nevts_, eventExtractor->n_events()); 

    for (int i=0;i<nevent;++i) 
    {
      if (i%10000 == 0)
        std::cout << "Processing " << i << "th event" << std::endl;

      PatExtractor::getInfo(i);// Retrieve the info from an existing ROOTuple      
      // Execute each plugins
      for (auto& plugin: m_plugins) {
        plugin->analyze(setup, *this);
      }

      ++nevent_tot; 
    }
  }
}


// What to do for each event

void PatExtractor::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  using namespace reco;

  if (do_fill_) 
  {
    PatExtractor::fillInfo(&event, setup); // Fill the ROOTuple
    // Execute each plugins
    for (auto& plugin: m_plugins) {
      plugin->analyze(event, setup, *this);
    }
  }

  ++nevent;
  ++nevent_tot;
}

void PatExtractor::endRun(Run const&, EventSetup const&) 
{
  std::cout << "Total # of events for this run   = "<< nevent  << std::endl;
}

void PatExtractor::endJob() {

  std::cout << "Total # of events for this job   = "<< nevent_tot     << std::endl;

  for (auto& extractor: m_extractors) {
    extractor->endJob(!do_fill_);
  }

  m_outfile->cd();

  for (auto& plugin: m_plugins) {
    plugin->endJob();
  }

  if (do_fill_) 
  {
    m_outfile->Write();
    m_outfile->Close();
  }
  else
  {
    m_infile->Close();
    m_outfile->Write();
    m_outfile->Close();
  }

}


// Here we fill the rootuple with info coming from the PatTuple

void PatExtractor::fillInfo(const edm::Event *event, const edm::EventSetup& iSetup) 
{
  MCExtractor* mcExtractor = static_cast<MCExtractor*>(getExtractor("MC").get());
  for (auto& extractor: m_extractors) {
    try {
      extractor->writeInfo(*event, iSetup, mcExtractor);
    } catch (cms::Exception& e) {
      std::stringstream context;
      context << "Calling event method for extractor " << typeid(*extractor).name() << "/\'" << extractor->getName() << "\'";
      e.addContext(context.str());
      throw e;
    }
  }
}


// Here we retrieve the info from an existing extracted ROOTuple 

void PatExtractor::getInfo(int ievent) 
{
  for (auto& extractor: m_extractors) {
    if (extractor->isOK())
      extractor->getInfo(ievent);
  }
}



// Here are the initializations when starting from scratch (need to create the extracted stuff)

void PatExtractor::initialize(const edm::ParameterSet& config) 
{
  m_outfile  = TFile::Open(outFilename_.c_str(),"RECREATE");
 
  // Load plugins
  if (!config.existsAs<edm::ParameterSet>("extractors")) {
    throw new std::logic_error("No extractors specified");
  }

  std::cout << std::endl << "Extractors: " << std::endl;
  const edm::ParameterSet& extractors = config.getParameterSet("extractors");
  std::vector<std::string> extractorNames = extractors.getParameterNames();
  for (std::string& extractorName: extractorNames) {
    edm::ParameterSet extractorData = extractors.getParameterSet(extractorName);
    bool enable = extractorData.getParameter<bool>("enable");
    if (! enable)
      continue;

    const std::string type = extractorData.getParameter<std::string>("type");
    edm::ParameterSet extractorParameters;
    if (extractorData.existsAs<edm::ParameterSet>("parameters"))
      extractorParameters = extractorData.getParameterSet("parameters");

    std::cout << " -> Adding extractor '" << extractorName << "' of type '" << type << "'" << std::endl;
    auto extractor = std::shared_ptr<SuperBaseExtractor>(PatExtractorExtractorFactory::get()->create(type, extractorName, extractorParameters));
    extractor->setIsMC(is_MC_);
    extractor->doConsumes(consumesCollector());
    extractor->check();

    addExtractor(extractorName, extractor);
  }
  std::cout << std::endl;
}

// Here are the initializations when starting from already extracted stuff

void PatExtractor::retrieve(const edm::ParameterSet& config) 
{
  m_infile     = TFile::Open(inFilename_.c_str(), "READ");
  m_outfile = TFile::Open(outFilename_.c_str(), "RECREATE");
  
  // Load plugins
  if (!config.existsAs<edm::ParameterSet>("extractors")) {
    throw new std::logic_error("No extractors specified");
  }

  std::cout << std::endl << "Extractors: " << std::endl;
  const edm::ParameterSet& extractors = config.getParameterSet("extractors");
  std::vector<std::string> extractorNames = extractors.getParameterNames();
  for (std::string& extractorName: extractorNames) {
    edm::ParameterSet extractorData = extractors.getParameterSet(extractorName);
    bool enable = extractorData.getParameter<bool>("enable");
    if (! enable)
      continue;

    const std::string type = extractorData.getParameter<std::string>("type");
    edm::ParameterSet extractorParameters;
    if (extractorData.existsAs<edm::ParameterSet>("parameters"))
      extractorParameters = extractorData.getParameterSet("parameters");

    std::cout << " -> Adding extractor '" << extractorName << "' of type '" << type << "'" << std::endl;
    auto extractor = std::shared_ptr<SuperBaseExtractor>(PatExtractorExtractorReadOnlyFactory::get()->create(type, extractorName, extractorParameters, m_infile));
    extractor->setIsMC(is_MC_);

    addExtractor(extractorName, extractor);
  }
}

DEFINE_FWK_MODULE(PatExtractor);
