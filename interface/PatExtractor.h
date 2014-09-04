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
#include "../interface/EventExtractor.h"
#include "../interface/ScaleFactorService.h"

#include <Extractors/PatExtractor/interface/ExtractorPlugin.h>

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
    if (m_extractorsIndexes.count(name))
      return m_extractors[m_extractorsIndexes[name]];
    else
      return std::shared_ptr<SuperBaseExtractor>();
  }

 private:

  void addExtractor(const std::string& name, const std::shared_ptr<SuperBaseExtractor>& extractor) {
    m_extractors.push_back(extractor);
    m_extractorsIndexes[name] = m_extractors.size() - 1;
  }

  bool is_MC_;
  bool do_fill_;
  int  nevts_;

  // Definition of root-tuple :

  std::string outFilename_;
  std::string inFilename_;

  TFile* m_infile;
  TFile* m_outfile;

  std::vector<std::shared_ptr<SuperBaseExtractor>> m_extractors;
  std::map<std::string, size_t> m_extractorsIndexes;

  std::vector<std::shared_ptr<patextractor::Plugin>> m_plugins;

  int iseventselected;
};


#endif
