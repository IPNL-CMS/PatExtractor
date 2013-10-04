#include "../interface/PatExtractor.h"

using namespace std;
using namespace edm;

PatExtractor::PatExtractor(const edm::ParameterSet& config) :
  is_MC_         (config.getUntrackedParameter<bool>("isMC", true)),
  do_fill_       (config.getUntrackedParameter<bool>("fillTree", true)),
  do_HLT_        (config.getUntrackedParameter<bool>("doHLT", false)),
  do_MC_         (config.getUntrackedParameter<bool>("doMC", false)),
  do_MCjpsi_     (config.getUntrackedParameter<bool>("doMCjpsi", false)),
  do_Photon_     (config.getUntrackedParameter<bool>("doPhoton", false)),
  do_Electron_   (config.getUntrackedParameter<bool>("doElectron", false)),

  do_Jet_        (config.getUntrackedParameter<bool>("doJet", false)),
  do_Muon_       (config.getUntrackedParameter<bool>("doMuon", false)),
  do_MET_        (config.getUntrackedParameter<bool>("doMET", false)),
  do_Vertex_     (config.getUntrackedParameter<bool>("doVertex", false)),
  do_Trk_        (config.getUntrackedParameter<bool>("doTrack", false)),
  do_PF_         (config.getUntrackedParameter<bool>("doPF", false)),
  do_Mtt_        (config.getUntrackedParameter<bool>("doMtt", false)),
  do_dimu_       (config.getUntrackedParameter<bool>("doDimuon", false)),
  do_ftt_        (config.getUntrackedParameter<bool>("do4TopHLT", false)),
  nevts_         (config.getUntrackedParameter<int>("n_events", 10000)),

  photon_tag_    (config.getParameter<edm::InputTag>("photon_tag")),
  electron_tag_  (config.getParameter<edm::InputTag>("electron_tag")),
  muon_tag_      (config.getParameter<edm::InputTag>("muon_tag")),
  MC_tag_        (config.getParameter<edm::InputTag>("MC_tag")),
  vtx_tag_       (config.getParameter<edm::InputTag>("vtx_tag")),
  trk_tag_       (config.getParameter<edm::InputTag>("trk_tag")),
  pf_tag_        (config.getParameter<edm::InputTag>("pf_tag")),

  outFilename_   (config.getParameter<std::string>("extractedRootFile")),
  inFilename_    (config.getParameter<std::string>("inputRootFile")),
  m_settings_    (config.getUntrackedParameter<std::vector<std::string> >("analysisSettings"))
{
  // We parse the analysis settings
  m_ana_settings = new AnalysisSettings(&m_settings_);
  m_ana_settings->parseSettings();

  if (do_Mtt_ && config.exists("mtt")) {
    m_mttParameterSet = config.getParameter<edm::ParameterSet>("mtt");
  }

  // Initializations

  nevent_tot = 0;

  // If do_fill is set to True, you extract the whole data, otherwise you start 
  // from a file already extracted (inFilename_)

  m_scaleFactors.reset(new ScaleFactorService(config));

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

  // Here is the small example analysis (dimuon mass spectra)

  if (do_dimu_ && do_Muon_)      
    m_dimuon_analysis = new dimuon_analysis(m_ana_settings);


  if (do_ftt_ && do_HLT_ && do_MC_)      
    m_fourtop_trigger_analysis = new fourtop_trigger_analysis(m_ana_settings);
}


// Job initialization
void PatExtractor::beginJob()
{
}


// What to do at the start of a new run
void PatExtractor::beginRun(Run const& run, EventSetup const& setup)
{

  if (is_MC_)
    m_scaleFactors->prepareBTaggingScaleFactors(setup);

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

  if (do_ftt_&& do_HLT_ && do_MC_) 
    m_fourtop_trigger_analysis->fourtop_trigger_finalize(nevent_tot);

}


// Here we fill the rootuple with info coming from the PatTuple

void PatExtractor::fillInfo(const edm::Event *event, const edm::EventSetup& iSetup) 
{
  // m_event->writeInfo(event,do_MC_);

  //if (do_HLT_)      m_HLT->writeInfo(event);
  //if (do_MET_)      m_MET->writeInfo(event);
  //if (do_Vertex_)   m_vertex->writeInfo(event);
  //if (do_Trk_)      m_track->writeInfo(event);
  //if (do_MC_)       m_MC->writeInfo(event);
  //if (do_Electron_) m_electron->writeInfo(event, m_MC, do_MC_);
  //if (do_Muon_)     m_muon->writeInfo(event, m_MC, do_MC_);
  //if (do_Jet_)      m_jet->writeInfo(event, iSetup, m_MC, do_MC_);

  //if (do_Photon_)   m_photon->writeInfo(event, m_MC, do_MC_);

  MCExtractor* mcExtractor = nullptr;
  if (do_MC_)
    mcExtractor = static_cast<MCExtractor*>(getExtractor("MC").get());

  for (auto& extractor: m_extractors)
    extractor->writeInfo(*event, iSetup, mcExtractor);
}


// Here we retrieve the info from an existing extracted ROOTuple 

void PatExtractor::getInfo(int ievent) 
{
  //m_event->getInfo(ievent);

  //if (do_HLT_)      m_HLT->getInfo(ievent);
  //if (do_MC_)       m_MC->getInfo(ievent);
  //if (do_Trk_)      m_track->getInfo(ievent);
  //if (do_Vertex_)   m_vertex->getInfo(ievent);
  //if (do_MET_)      m_MET->getInfo(ievent);
  //if (do_Muon_)     m_muon->getInfo(ievent);
  //if (do_Electron_) m_electron->getInfo(ievent);
  //if (do_Jet_)      m_jet->getInfo(ievent);
  //if (do_Photon_)   m_photon->getInfo(ievent);

  for (auto& extractor: m_extractors) {
    if (extractor->isOK())
      extractor->getInfo(ievent);
  }
}



// Here are the initializations when starting from scratch (need to create the extracted stuff)

void PatExtractor::initialize(const edm::ParameterSet& config) 
{
  m_outfile  = TFile::Open(outFilename_.c_str(),"RECREATE");
  //m_event    = new EventExtractor();

  //m_HLT      = new HLTExtractor(do_HLT_);
  //m_MC       = new MCExtractor(do_MC_);
  //m_photon   = new PhotonExtractor(do_Photon_,photon_tag_);
  //m_electron = new ElectronExtractor(do_Electron_,electron_tag_);
  //m_jet      = new JetExtractor(do_Jet_, jet_tag_, correctJets_, jetCorrectorLabel_);
  //m_MET      = new METExtractor(do_MET_,met_tag_);
  //m_muon     = new MuonExtractor(do_Muon_,muon_tag_);
  //m_vertex   = new VertexExtractor(do_Vertex_,vtx_tag_);
  //m_track    = new TrackExtractor(do_Trk_,trk_tag_);

  // Register extractors
  addExtractor("event", new EventExtractor("event"));

  if (do_MC_)
    addExtractor("MC", new MCExtractor("MC", do_MC_, do_MCjpsi_));

  if (do_HLT_)
    addExtractor("HLT", new HLTExtractor("HLT", do_HLT_, config));

  if (do_Trk_)
    addExtractor("track", new TrackExtractor("track", trk_tag_, do_Trk_));

  if (do_PF_)
    addExtractor("PFpart", new PFpartExtractor("PFpart", pf_tag_, do_PF_));

  if (do_Vertex_)
    addExtractor("vertex", new VertexExtractor("Vertices", vtx_tag_, do_Vertex_));

  if (do_Electron_)
    addExtractor("electrons", new ElectronExtractor("electron_PF", electron_tag_, do_Electron_));

  if (do_Muon_)
    addExtractor("muons", new MuonExtractor("muon_PF", muon_tag_, vtx_tag_, do_Muon_));

  if (do_Jet_ || do_MET_)
    addExtractor("JetMET", new JetMETExtractor("jet_PF", "MET_PF", config));

  if (do_Photon_)
    addExtractor("photons", new PhotonExtractor("photon", photon_tag_, do_Photon_));

  // Add non isolated leptons for vetoes
  if (do_Electron_)
    addExtractor("electrons_loose", new ElectronExtractor("electron_loose_PF", edm::InputTag("selectedPatElectronsLoosePFlow"), true));

  if (do_Muon_)
    addExtractor("muons_loose", new MuonExtractor("muon_loose_PF", edm::InputTag("selectedPatMuonsLoosePFlow"), vtx_tag_, true));

  for (auto& extractor: m_extractors) {
    extractor->setIsMC(is_MC_);
    extractor->setScaleFactorsService(m_scaleFactors);
  }
}




// Here are the initializations when starting from already extracted stuff

void PatExtractor::retrieve(const edm::ParameterSet& config) 
{
  m_infile     = TFile::Open(inFilename_.c_str(), "READ");
  m_outfile    = TFile::Open(outFilename_.c_str(), "RECREATE");

  // AOD content
  //m_event      = new EventExtractor(m_infile);
  //m_HLT        = new HLTExtractor(m_infile);
  //m_MC         = new MCExtractor(m_infile);
  //m_vertex     = new VertexExtractor(m_infile);
  //m_track      = new TrackExtractor(m_infile);

  // PAT content
  //m_MET        = new METExtractor(m_infile);
  //m_muon       = new MuonExtractor(m_infile);
  //m_photon     = new PhotonExtractor(m_infile);
  //m_electron   = new ElectronExtractor(m_infile);
  //m_jet        = new JetExtractor(m_infile);

  // Register extractors
  addExtractor("event", new EventExtractor("event", m_infile));
  addExtractor("MC", new MCExtractor("MC", m_infile));
  addExtractor("HLT", new HLTExtractor("HLT", m_infile));
  addExtractor("track", new TrackExtractor("track", m_infile));
  addExtractor("PFpart", new PFpartExtractor("PFpart", m_infile));

  addExtractor("vertex", new VertexExtractor("Vertices", m_infile));

  addExtractor("electrons", new ElectronExtractor("electron_PF", m_infile));
  addExtractor("electrons_loose", new ElectronExtractor("electron_loose_PF", m_infile));

  addExtractor("muons", new MuonExtractor("muon_PF", m_infile));
  addExtractor("muons_loose", new MuonExtractor("muon_loose_PF", m_infile));

  addExtractor("JetMET", new JetMETExtractor("jet_PF", "MET_PF", m_infile));
  addExtractor("photons", new PhotonExtractor("photon", m_infile));

  // We set some variables wrt the info retrieved (if the tree is not there, don't go further...)  
  do_HLT_      = getExtractor("HLT")->isOK();
  do_MC_       = getExtractor("MC")->isOK();
  do_Photon_   = getExtractor("photons")->isOK();
  do_Electron_ = getExtractor("electrons")->isOK();
  do_Jet_      = getExtractor("JetMET")->isOK();
  do_Muon_     = getExtractor("muons")->isOK();
  do_Vertex_   = getExtractor("vertex")->isOK();
  do_Trk_      = getExtractor("track")->isOK();
  do_PF_       = getExtractor("PFpart")->isOK();

  for (auto& extractor: m_extractors) {
    extractor->setIsMC(is_MC_);
    extractor->setScaleFactorsService(m_scaleFactors);
  }
}
