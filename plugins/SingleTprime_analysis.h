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
    int patIndexToExtractorIndex(int patIndex) const;

    void MCidentification();
    void reset();
    void fillTree();
    
  private:
    int j;

    TTree* m_tree_stp;
    
    int m_evt;
    int m_nPU; //PileUp
    int evt_num;
    
    float m_THT;
    float m_jet1pt;
    float m_jet2pt;
    float m_jet3pt;
    float m_jet4pt;
    float m_jet5pt;
    float m_jet6pt;
    float m_DRHiggsJets;
    float m_DRWJets;
    float m_DRTopHiggs;
    float m_DRWHiggs;
    float m_RelTHT;
    float m_Sphericity;
    float m_Aplanarity;
    float m_DRTrueHiggsRecoHiggs;
    float m_DRTrueWRecoW;
    float m_DRTrueTopRecoTop;
    float m_DRTrueTprimeRecoTprime;
    float m_DRTrueFirstHiggsJetRecoJet;
    float m_DRTrueSecondHiggsJetRecoJet;
    float m_DRTrueFirstWJetRecoJet;
    float m_DRTrueSecondWJetRecoJet; 
    float m_DRTrueTopJetRecoJet;   

    float m_DRTrueWJets;
    float m_DRMatchedWJets;
    float m_DPhiTrueWJets;   
    float m_DPhiMatchedWJets;

    //bool isMC;
    float m_jet_Ptcut;
    float m_jet_EtaMaxcut;
    float m_jet_EtaAccepcut;
    float m_jet_OverlapAccep;
    float m_jet_MultInAcceptance;
    float m_jet_MultOutAcceptance;
    float m_DRMatching;
    float m_DPtMatching;

    //Booleans to activate each cut
    bool Cut0;
    bool Cut1;
    bool Cut2;
    bool Cut3;
    bool Cut4;
    bool Cut5;
    bool Cut6;
    bool Cut7;
    bool Cut8;
    bool Cut9;
    bool Cut10;
    bool Cut11;
    bool Cut12;
    bool Cut13;

    bool DoMCMatching;
 
    //Reconstructed particles
    TLorentzVector* ReconstructedHiggs;
    TLorentzVector* ReconstructedW;
    TLorentzVector* ReconstructedTop;
    TLorentzVector* ReconstructedTprime;
    TLorentzVector* TrueHiggs;
    TLorentzVector* TrueW;
    TLorentzVector* TrueTop;
    TLorentzVector* TrueTprime;
    TLorentzVector* TrueTprimeAcompainingJet;
    TLorentzVector* FirstHiggsJet;
    TLorentzVector* SecondHiggsJet;
    TLorentzVector* FirstWJet;
    TLorentzVector* SecondWJet;
    TLorentzVector* TopJet;
    TLorentzVector* FirstTrueHiggsJet;
    TLorentzVector* SecondTrueHiggsJet;
    TLorentzVector* FirstTrueWJet;
    TLorentzVector* SecondTrueWJet;
    TLorentzVector* TopTrueJet;

    //Matching with MC truth
    int CorrectTprime;
    int CorrectH;
    int CorrectW;
    int CorrectTop;
    int CorrectHiggsJet;
    int CorrectWJet;
    int CorrectTopJet;
    int NumbMatchedHiggsJets;
    int NumbMatchedWJets;
    int NumbMatchedTopJets;

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

    /// scale factors for jet energy resolution
    std::vector<double> jetEnergyResolutionScaleFactors_;
    std::vector<double> jetEnergyResolutionEtaBinning_;

    bool m_trigger_passed;	

    float m_weight;	
    float m_weight_error_low;
    float m_weight_error_high;

  };

}
