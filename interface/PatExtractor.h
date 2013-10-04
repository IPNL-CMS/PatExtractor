#ifndef PatExtractor_h
#define PatExtractor_h

/** \class PatExtractor
 *  Class that produces a roottuple from PATuples
 */

#include <memory>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "../interface/BaseExtractor.h"
#include "../interface/ElectronExtractor.h"
#include "../interface/MCExtractor.h"
#include "../interface/PhotonExtractor.h"
#include "../interface/JetMETExtractor.h"
#include "../interface/METExtractor.h"
#include "../interface/MuonExtractor.h"
#include "../interface/VertexExtractor.h"
#include "../interface/EventExtractor.h"
#include "../interface/HLTExtractor.h"
#include "../interface/AnalysisSettings.h"
#include "../interface/TrackExtractor.h"
#include "../interface/PFpartExtractor.h"
#include "../interface/ScaleFactorService.h"

#include <Extractors/PatExtractor/interface/ExtractorPlugin.h>

#include "../interface/dimuon_analysis.h"
#include "../interface/fourtop_trigger_analysis.h"

#include "TFile.h"

class PatExtractor : public edm::EDAnalyzer {
 public:
  /// Constructor
  PatExtractor(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~PatExtractor(){ }
  int nevent;
  int nevent_tot;

  virtual void beginJob();
  virtual void endJob();

  void beginRun(edm::Run const&, edm::EventSetup const&);
  void endRun(edm::Run const&, edm::EventSetup const&);

  void analyze(const edm::Event&, const edm::EventSetup& );
  
  void fillInfo(const edm::Event *event, const edm::EventSetup& iSetup);
  void getInfo(int ievent);
  void initialize(const edm::ParameterSet&);
  void retrieve(const edm::ParameterSet&);

  std::shared_ptr<SuperBaseExtractor> getExtractor(const std::string& name) {
    return m_extractors[m_extractorsIndexes[name]];
  }

 private:

  void addExtractor(const std::string& name, SuperBaseExtractor* extractor) {
    m_extractors.push_back(std::shared_ptr<SuperBaseExtractor>(extractor));
    m_extractorsIndexes[name] = m_extractors.size() - 1;
  }

  bool is_MC_;

  bool do_fill_;
  bool do_HLT_;
  bool do_MC_;
  bool do_MCjpsi_;
  bool do_Photon_;
  bool do_Electron_;

  bool do_Jet_;
  bool do_Muon_;
  bool do_MET_;
  bool do_Vertex_;
  bool do_Trk_;
  bool do_PF_;

  bool do_Mtt_;
  bool do_dimu_;
  bool do_ftt_;
  int  nevts_;

  edm::InputTag photon_tag_;   // 
  edm::InputTag electron_tag_; // 
  edm::InputTag muon_tag_;     // 
  edm::InputTag MC_tag_;       // 
  edm::InputTag vtx_tag_;      // 
  edm::InputTag trk_tag_;      // 
  edm::InputTag pf_tag_;      //

  // Definition of root-tuple :

  std::string outFilename_;
  std::string inFilename_;

  std::vector<std::string> m_settings_;

  edm::ParameterSet m_mttParameterSet;

  TFile* m_infile;
  TFile* m_outfile;

  std::vector<std::shared_ptr<SuperBaseExtractor>> m_extractors;
  std::map<std::string, size_t> m_extractorsIndexes;

  std::vector<std::shared_ptr<patextractor::Plugin>> m_plugins;

  AnalysisSettings*  m_ana_settings;

  dimuon_analysis*            m_dimuon_analysis;
  fourtop_trigger_analysis*   m_fourtop_trigger_analysis;

  int iseventselected;

  std::shared_ptr<ScaleFactorService> m_scaleFactors;
};


#endif
