#include <Extractors/PatExtractor/interface/ExtractorPlugin.h>

#include <vector>
#include <boost/regex.hpp>

#include <TTree.h>
#include <TLorentzVector.h>
#include <TRef.h>

class AnalysisSettings;
class EventExtractor;
class MuonExtractor;
class ElectronExtractor;
class METExtractor;
class VertexExtractor;
class KinFit;
class HLTExtractor;
class PatExtractor;
class MCExtractor;
class JetMETExtractor;

namespace edm {
  class EventSetup;
}

namespace patextractor {

  //class SingleTprime_analysis: patextractor::Plugin {
  class SingleTprime_analysis: public Plugin {
  public:
    SingleTprime_analysis(const edm::ParameterSet& iConfig);
    ~SingleTprime_analysis();
    
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, PatExtractor& extractor);
    virtual void analyze(const edm::EventSetup& iSetup, PatExtractor& extractor);

    //Functions
    bool isJetSel(TLorentzVector *jet);
    bool isJetAccepSel(TLorentzVector *jet);
    bool isJetForwSel(TLorentzVector *jet);
    int  SingleTprime_Sel();

  void reset();
  void fillTree();
    
  private:
    int j;

    TTree* m_tree_stp;
    
    int m_evt;
    
    float m_THT;
    float m_jet1pt;
    float m_jet2pt;
    float m_jet3pt;
    float m_jet4pt;
    float m_jet5pt;
    float m_jet6pt;
    
    //bool isMC;
    float m_jet_Ptcut;
    float m_jet_EtaMaxcut;
    float m_jet_EtaAccepcut;
    float m_jet_MultInAcceptance;
    float m_jet_MultOutAcceptance;

    //Linking extractors
    std::shared_ptr<EventExtractor> m_event;
    std::shared_ptr<MCExtractor>    m_MC;
    std::shared_ptr<VertexExtractor> m_vertex;
    std::shared_ptr<JetMETExtractor> m_jetMet;
    float m_JET_btag_CSVL;
    float m_JET_btag_CSVM;
    float m_JET_btag_CSVT;

    //Counter for b-tagged jets
    int m_NBtaggedJets_CSVL;

  };

}
