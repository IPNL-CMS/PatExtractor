#include <Extractors/PatExtractor/interface/ExtractorPlugin.h>

#include <TH2F.h>
#include <vector>
#include <boost/regex.hpp>

#include <Extractors/PatExtractor/interface/ExtractorPlugin.h>

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

class SingleTprime: patextractor::Plugin {
 public:
  SingleTprime(const edm::ParameterSet& iConfig);
  
  virtual void analyze(const edm::EventSetup& iSetup, PatExtractor& extractor);

 private:
  TH1F *M12345;
  TH1F *Cut0;
  TH1F *Cut1;
  TH1F *Cut2;
  TH1F *Cut3;
  TH1F *Cut4;
  TH1F *Cut5;
  TH1F *Cut6;
  TH1F *Cut7;
  TH1F *Cut8;
  TH1F *Cut9;
  TH1F *Cut10;
  TH1F *Cut11;
  TH1F *Cut12;
  TH1F *Cut13;
  TH1F *Cut14;
  TH1F *Jet1_pt;
  TH1F *Jet2_pt;
  TH1F *Jet3_pt;
  TH1F *Jet4_pt;
  TH1F *Jet5_pt;
  TH1F *Jet6_pt;
  TH1F *THT;
  TH1F *DeltaRH;
  TH1F *MAroundHiggs;
  TH1F *MAroundW;
  TH1F *MAroundTop;
  TH1F *jetnum;
  TH1F *HiggsJetsMass;
  TH1F *RelHT;
  TH1F *AplanarityHad;
  TH2F *MHDR;
  TH2F *HptToppt;
  TH2F *DRWDRH;
  TH2F *DPHDPT;
  TH2F *DPHDPW;
  TH2F *TprimeMassVsHiggsMass;
  TH2F *TprimeMassVsDRTopHiggs;
  TH2F *TprimeMassVsTprimePT;
  
};
