#include "../interface/ElectronExtractor.h"

#include <DataFormats/VertexReco/interface/Vertex.h>
#include <EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h>
#include <EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h>
#include <DataFormats/RecoCandidate/interface/IsoDeposit.h>
#include <DataFormats/RecoCandidate/interface/IsoDepositVetos.h>

ElectronExtractor::ElectronExtractor(const std::string& name, std::shared_ptr<ScaleFactorService> sf, const edm::InputTag& tag, bool doTree)
  : BaseExtractor(name, sf), m_ele_lorentzvector(nullptr)
{
  m_tag = tag;

  m_ele_lorentzvector = new TClonesArray("TLorentzVector");

  const auto& sfWorkingPoints = m_scaleFactorService->getElectronScaleFactorWorkingPoints();
  for (auto& it: sfWorkingPoints) {
      std::string name = "electron_scaleFactor_" + ScaleFactorService::workingPointToString(it.first) + "eff_" + ScaleFactorService::workingPointToString(it.second) + "iso";
      m_scaleFactors[name] = ScaleFactorCollection();
      m_scaleFactors[name].setWriteMode();
  };

  setPF((tag.label()).find("PFlow"));
  reset();

  // Tree definition

  m_OK = false;

  if (doTree) {
    m_OK = true;

    m_tree_electron     = new TTree(m_name.c_str(), "PAT PF electron info");  
    m_tree_electron->Branch("n_electrons",                       &m_size, "n_electrons/i");  
    m_tree_electron->Branch("electron_4vector","TClonesArray",   &m_ele_lorentzvector, 1000, 0);
    m_tree_electron->Branch("electron_vx",                       &m_ele_vx,     "electron_vx[n_electrons]/F");  
    m_tree_electron->Branch("electron_vy",                       &m_ele_vy,     "electron_vy[n_electrons]/F");  
    m_tree_electron->Branch("electron_vz",                       &m_ele_vz,     "electron_vz[n_electrons]/F");  
    m_tree_electron->Branch("electron_charge",                   &m_ele_charge,    "electron_charge[n_electrons]/I");  

    m_tree_electron->Branch("electron_dB",                       &m_ele_dB,        "electron_dB[n_electrons]/F");  
    m_tree_electron->Branch("electron_trackIso",                 &m_ele_trackIso,  "electron_trackIso[n_electrons]/F");  
    m_tree_electron->Branch("electron_ecalIso",                  &m_ele_ecalIso,   "electron_ecalIso[n_electrons]/F");  
    m_tree_electron->Branch("electron_hcalIso",                  &m_ele_hcalIso,   "electron_hcalIso[n_electrons]/F");  

    m_tree_electron->Branch("electron_pfParticleIso",            &m_ele_pfParticleIso,     "electron_pfParticleIso[n_electrons]/F");
    m_tree_electron->Branch("electron_pfChargedHadronIso",       &m_ele_pfChargedHadronIso,"electron_pfChargedHadronIso[n_electrons]/F");
    m_tree_electron->Branch("electron_pfNeutralHadronIso",       &m_ele_pfNeutralHadronIso,"electron_pfNeutralHadronIso[n_electrons]/F");
    m_tree_electron->Branch("electron_pfPhotonIso",              &m_ele_pfPhotonIso,       "electron_pfPhotonIso[n_electrons]/F");

    m_tree_electron->Branch("electron_numberOfMissedInnerLayer", &m_ele_numberOfMissedInnerLayer, "electron_numberOfMissedInnerLayer[n_electrons]/I");  
    m_tree_electron->Branch("electron_mcParticleIndex",          &m_ele_MCIndex,   "electron_mcParticleIndex[n_electrons]/I");  

    m_tree_electron->Branch("electron_eidMVATrigV0",             &m_ele_eidMVATrigV0, "electron_eidMVATrigV0[n_electrons]/F");
    m_tree_electron->Branch("electron_passConversionVeto",       &m_ele_passConversionVeto, "electron_passConversionVeto[n_electrons]/O");
    m_tree_electron->Branch("electron_SCEta",                    &m_ele_SCEta, "electron_SCEta[n_electrons]/F");
    m_tree_electron->Branch("electron_effectiveArea",            &m_ele_effectiveArea, "electron_effectiveArea[n_electrons]/F");
    m_tree_electron->Branch("electron_relIsolation",             &m_ele_relIsolation, "electron_relIsolation[n_electrons]/F");
    m_tree_electron->Branch("electron_rhoCorrectedRelIsolation", &m_ele_rhoCorrectedRelIsolation, "electron_rhoCorrectedRelIsolation[n_electrons]/F");
    m_tree_electron->Branch("electron_deltaBetaCorrectedRelIsolation", &m_ele_deltaBetaCorrectedRelIsolation, "electron_deltaBetaCorrectedRelIsolation[n_electrons]/F");

    m_tree_electron->Branch("electron_passVetoID",               &m_ele_passVetoID, "electron_passVetoID[n_electrons]/O");
    m_tree_electron->Branch("electron_passLooseID",              &m_ele_passLooseID, "electron_passLooseID[n_electrons]/O");
    m_tree_electron->Branch("electron_passMediumID",             &m_ele_passMediumID, "electron_passMediumID[n_electrons]/O");
    m_tree_electron->Branch("electron_passTightID",              &m_ele_passTightID, "electron_passTightID[n_electrons]/O");

    for (auto& it: m_scaleFactors) {
      m_tree_electron->Branch(it.first.c_str(), & it.second.getBackingArray());
    }
  }
}

  ElectronExtractor::ElectronExtractor(const std::string& name, std::shared_ptr<ScaleFactorService> sf, TFile* a_file)
:  BaseExtractor(name, sf), m_ele_lorentzvector(nullptr)
{
  m_file = a_file;

  // Tree definition
  m_OK = false;

  m_tree_electron = dynamic_cast<TTree*>(a_file->Get(m_name.c_str()));

  if (!m_tree_electron)
  {
    std::cerr << "This tree doesn't exist!!!" << std::endl;
    return;  
  }

  m_OK = true;

  m_ele_lorentzvector = new TClonesArray("TLorentzVector");

  // Branches definition (with safety check)

  if (m_tree_electron->FindBranch("n_electrons")) 
    m_tree_electron->SetBranchAddress("n_electrons",&m_size);
  if (m_tree_electron->FindBranch("electron_4vector"))
    m_tree_electron->SetBranchAddress("electron_4vector",&m_ele_lorentzvector);
  if (m_tree_electron->FindBranch("electron_vx"))
    m_tree_electron->SetBranchAddress("electron_vx",&m_ele_vx);
  if (m_tree_electron->FindBranch("electron_vy"))
    m_tree_electron->SetBranchAddress("electron_vy",&m_ele_vy);
  if (m_tree_electron->FindBranch("electron_vz"))
    m_tree_electron->SetBranchAddress("electron_vz",&m_ele_vz);
  if (m_tree_electron->FindBranch("electron_charge"))
    m_tree_electron->SetBranchAddress("electron_charge",&m_ele_charge);

  if (m_tree_electron->FindBranch("electron_dB"))
    m_tree_electron->SetBranchAddress("electron_dB",&m_ele_dB);
  if (m_tree_electron->FindBranch("electron_trackIso"))
    m_tree_electron->SetBranchAddress("electron_trackIso",&m_ele_trackIso);
  if (m_tree_electron->FindBranch("electron_ecalIso"))
    m_tree_electron->SetBranchAddress("electron_ecalIso",&m_ele_ecalIso);
  if (m_tree_electron->FindBranch("electron_hcalIso"))
    m_tree_electron->SetBranchAddress("electron_hcalIso",&m_ele_hcalIso);
  if (m_tree_electron->FindBranch("electron_pfParticleIso"))
    m_tree_electron->SetBranchAddress("electron_pfParticleIso",&m_ele_pfParticleIso);
  if (m_tree_electron->FindBranch("electron_pfChargedHadronIso"))
    m_tree_electron->SetBranchAddress("electron_pfChargedHadronIso",&m_ele_pfChargedHadronIso);
  if (m_tree_electron->FindBranch("electron_pfNeutralHadronIso"))
    m_tree_electron->SetBranchAddress("electron_pfNeutralHadronIso",&m_ele_pfNeutralHadronIso);
  if (m_tree_electron->FindBranch("electron_pfPhotonIso"))
    m_tree_electron->SetBranchAddress("electron_pfPhotonIso",&m_ele_pfPhotonIso);
  if (m_tree_electron->FindBranch("electron_numberOfMissedInnerLayer"))
    m_tree_electron->SetBranchAddress("electron_numberOfMissedInnerLayer",&m_ele_numberOfMissedInnerLayer);
  if (m_tree_electron->FindBranch("electron_mcParticleIndex"))
    m_tree_electron->SetBranchAddress("electron_mcParticleIndex",&m_ele_MCIndex);

  if (m_tree_electron->FindBranch("electron_eidMVATrigV0"))
    m_tree_electron->SetBranchAddress("electron_eidMVATrigV0", &m_ele_eidMVATrigV0);

  if (m_tree_electron->FindBranch("electron_passConversionVeto"))
    m_tree_electron->SetBranchAddress("electron_passConversionVeto", &m_ele_passConversionVeto);

  if (m_tree_electron->FindBranch("electron_SCEta"))
    m_tree_electron->SetBranchAddress("electron_SCEta", &m_ele_SCEta);

  if (m_tree_electron->FindBranch("electron_passVetoID"))
    m_tree_electron->SetBranchAddress("electron_passVetoID", &m_ele_passVetoID);

  if (m_tree_electron->FindBranch("electron_passLooseID"))
    m_tree_electron->SetBranchAddress("electron_passLooseID", &m_ele_passLooseID);

  if (m_tree_electron->FindBranch("electron_passMediumID"))
    m_tree_electron->SetBranchAddress("electron_passMediumID", &m_ele_passMediumID);

  if (m_tree_electron->FindBranch("electron_passTightID"))
    m_tree_electron->SetBranchAddress("electron_passTightID", &m_ele_passTightID);

  if (m_tree_electron->FindBranch("electron_effectiveArea"))
    m_tree_electron->SetBranchAddress("electron_effectiveArea", &m_ele_effectiveArea);

  if (m_tree_electron->FindBranch("electron_relIsolation"))
    m_tree_electron->SetBranchAddress("electron_relIsolation", &m_ele_relIsolation);

  if (m_tree_electron->FindBranch("electron_rhoCorrectedRelIsolation"))
    m_tree_electron->SetBranchAddress("electron_rhoCorrectedRelIsolation", &m_ele_rhoCorrectedRelIsolation);

  if (m_tree_electron->FindBranch("electron_deltaBetaCorrectedRelIsolation"))
    m_tree_electron->SetBranchAddress("electron_deltaBetaCorrectedRelIsolation", &m_ele_deltaBetaCorrectedRelIsolation);

  for (auto& it: m_scaleFactors) {
    if (m_tree_electron->FindBranch(it.first.c_str()))
      m_tree_electron->Branch(it.first.c_str(), & it.second.getBackingArray());
  }
}


ElectronExtractor::~ElectronExtractor() {
  delete m_ele_lorentzvector;
}

void ElectronExtractor::beginJob(edm::ConsumesCollector&& collector, bool isInAnalysisMode) {
  BaseExtractor::beginJob(std::forward<edm::ConsumesCollector>(collector), isInAnalysisMode);

  m_allConversionsToken = collector.consumes<reco::ConversionCollection>(edm::InputTag("allConversions"));
  m_offlineBeamSpotToken = collector.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  m_primaryVerticesToken = collector.consumes<reco::VertexCollection>(edm::InputTag("goodOfflinePrimatryVertices"));

  edm::InputTag rhoInputTag;
  if (! m_isMC)
    rhoInputTag = edm::InputTag("kt6PFJets", "rho", "RECO");
  else
    rhoInputTag = edm::InputTag("kt6PFJetsForIsolation", "rho", "PAT");

  m_rhoToken = collector.consumes<double>(rhoInputTag);
}

//
// Method filling the main particle tree
//

void ElectronExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Electron& object, int index) {

  if (index >= m_electrons_MAX)
    return;

  const pat::Electron& part = static_cast<const pat::Electron&>(object);

  new((*m_ele_lorentzvector)[index]) TLorentzVector(part.px(),part.py(),part.pz(),part.energy());
  m_ele_vx[index]                 = part.vx();
  m_ele_vy[index]                 = part.vy();
  m_ele_vz[index]                 = part.vz();
  m_ele_charge[index]             = part.charge();
  m_ele_SCEta[index]              = part.superCluster()->eta();
  m_ele_passConversionVeto[index] = part.passConversionVeto();

  // 2012 Electron ID based on MVA
  if (part.isElectronIDAvailable("mvaTrigV0"))
    m_ele_eidMVATrigV0[index]        = part.electronID("mvaTrigV0");

  m_ele_dB[index]               = part.dB() ;
  m_ele_trackIso[index]         = part.trackIso();
  m_ele_ecalIso[index]          = part.ecalIso();
  m_ele_hcalIso[index]          = part.hcalIso();


  if (part.gsfTrack().isNonnull())
  {
    m_ele_numberOfMissedInnerLayer[index] = part.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
  }

  if (m_isPF)
  {
    m_ele_pfParticleIso[index]      = part.particleIso();
    m_ele_pfChargedHadronIso[index] = part.chargedHadronIso();
    m_ele_pfNeutralHadronIso[index] = part.neutralHadronIso();
    m_ele_pfPhotonIso[index]        = part.photonIso();
  }

  // See https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel#Electrons

  {
    // Check if this electron pass VETO criteria
    edm::Handle<reco::ConversionCollection> hConversions;
    event.getByToken(m_allConversionsToken, hConversions);

    // beam spot
    edm::Handle<reco::BeamSpot> hBeamspot;
    event.getByToken(m_offlineBeamSpotToken, hBeamspot);
    const reco::BeamSpot &beamSpot = *hBeamspot;

    // vertices
    edm::Handle<reco::VertexCollection> hVtx;
    event.getByToken(m_primaryVerticesToken, hVtx);

    edm::Handle<double> hRhoIso;
    event.getByToken(m_rhoToken, hRhoIso);
    double rhoIso = *hRhoIso;

    // Compute isolation in a cone of 0.3. One can use PAT functions chargedHadronIso(), neutralHadronIso() and photonIso(). They are supposed to do the same thing that what follow.
    // This is just to check that everything is consistant.
    //reco::isodeposit::AbsVetos chargedHadronVetos;
    //reco::isodeposit::AbsVetos photonVetos;

    //// PF isolation veto setup EGM recommendation
    //if (fabs(part.superCluster()->eta()) > 1.479 ) {
      //reco::isodeposit::Direction Dir = reco::isodeposit::Direction(part.superCluster()->eta(), part.superCluster()->phi());
      //chargedHadronVetos.push_back(new reco::isodeposit::ConeVeto(Dir, 0.015));
      //photonVetos.push_back(new reco::isodeposit::ConeVeto(Dir, 0.08));
    //}

    ////cone size 0.3 
    //const double chIso03 = part.isoDeposit(pat::PfChargedHadronIso)->depositWithin(0.3, chargedHadronVetos);
    //const double nhIso03 = part.isoDeposit(pat::PfNeutralHadronIso)->depositWithin(0.3);
    //const double phIso03 = part.isoDeposit(pat::PfGammaIso)->depositWithin(0.3, photonVetos);
    //const double puChIso03 = part.isoDeposit(pat::PfPUChargedHadronIso)->depositWithin(0.3, chargedHadronVetos);

    // All of the above can be replaced with
    const double chIso03 = part.chargedHadronIso();
    const double nhIso03 = part.neutralHadronIso();
    const double phIso03 = part.photonIso();
    const double puChIso03 = part.puChargedHadronIso();

    ElectronEffectiveArea::ElectronEffectiveAreaTarget ea = ElectronEffectiveArea::kEleEAData2012;
    if (m_isMC)
      ea = ElectronEffectiveArea::kEleEAFall11MC;

    m_ele_passVetoID[index] = EgammaCutBasedEleId::PassWP(
        EgammaCutBasedEleId::VETO,
        part,
        hConversions,
        beamSpot,
        hVtx,
        chIso03,
        phIso03,
        nhIso03,
        rhoIso);

    m_ele_passLooseID[index] = EgammaCutBasedEleId::PassWP(
        EgammaCutBasedEleId::LOOSE,
        part,
        hConversions,
        beamSpot,
        hVtx,
        chIso03,
        phIso03,
        nhIso03,
        rhoIso);

    m_ele_passMediumID[index] = EgammaCutBasedEleId::PassWP(
        EgammaCutBasedEleId::MEDIUM,
        part,
        hConversions,
        beamSpot,
        hVtx,
        chIso03,
        phIso03,
        nhIso03,
        rhoIso);

    m_ele_passTightID[index] = EgammaCutBasedEleId::PassWP(
        EgammaCutBasedEleId::TIGHT,
        part,
        hConversions,
        beamSpot,
        hVtx,
        chIso03,
        phIso03,
        nhIso03,
        rhoIso);

    // Compute effective area
    float AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, part.superCluster()->eta(), ea);

    m_ele_effectiveArea[index] = AEff03;

    const double relIso = (chIso03 + nhIso03 + phIso03) / part.pt();
    m_ele_relIsolation[index] = relIso;

    const double relIsoRhoCorrected = (chIso03 + std::max(0.0, (nhIso03 + phIso03) - rhoIso * AEff03)) / part.pt();
    m_ele_rhoCorrectedRelIsolation[index] = relIsoRhoCorrected;

    const double relIsoDeltaBetacCorrected = (chIso03 + std::max(0.0, (nhIso03 + phIso03) - 0.5 * puChIso03)) / part.pt();
    m_ele_deltaBetaCorrectedRelIsolation[index] = relIsoDeltaBetacCorrected;
  }

  if (m_isMC) {
    for (auto& it: m_scaleFactors) {
      std::pair<ScaleFactorService::WorkingPoint, ScaleFactorService::WorkingPoint> workingPoints = ScaleFactorService::getWorkingPointFromName(it.first);
      it.second.push_back(m_scaleFactorService->getElectronScaleFactor(workingPoints.first, workingPoints.second, part.pt(), part.eta()));
    }
  }
}

//
// Method getting the info from an input file
//

void ElectronExtractor::getInfo(int ievt) 
{
  m_tree_electron->GetEntry(ievt); 
}

// Method initializing everything (to do for each event)

void ElectronExtractor::reset()
{
  m_size = 0;
  for (auto& it: m_scaleFactors) {
    it.second.clear();
  }

  for (int i=0;i<m_electrons_MAX;++i) 
  {
    m_ele_vx[i] = 0.;
    m_ele_vy[i] = 0.;
    m_ele_vz[i] = 0.;
    m_ele_charge[i] = 0;

    m_ele_dB[i]=0;                          
    m_ele_trackIso[i]=0;                      
    m_ele_ecalIso[i]=0;
    m_ele_hcalIso[i]=0;
    m_ele_pfParticleIso[i]=0;
    m_ele_pfChargedHadronIso[i]=0;
    m_ele_pfNeutralHadronIso[i]=0;
    m_ele_pfPhotonIso[i]=0;
    m_ele_numberOfMissedInnerLayer[i]=0;
    m_ele_MCIndex[i]=-1;

    m_ele_eidMVATrigV0[i] = -1;
    m_ele_SCEta[i] = -1;
    m_ele_passConversionVeto[i] = false;
    m_ele_effectiveArea[i] = 0;
    m_ele_rhoCorrectedRelIsolation[i] = -1;
    m_ele_deltaBetaCorrectedRelIsolation[i] = -1;
    m_ele_relIsolation[i] = -1;
    m_ele_passVetoID[i] = false;
    m_ele_passLooseID[i] = false;
    m_ele_passMediumID[i] = false;
    m_ele_passTightID[i] = false;
  }

  m_ele_lorentzvector->Clear();
}


void ElectronExtractor::fillTree()
{
  m_tree_electron->Fill(); 
}
