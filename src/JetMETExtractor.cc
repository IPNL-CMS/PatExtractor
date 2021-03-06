#include "../interface/JetMETExtractor.h"
#include "../interface/JECReader.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/Math/interface/deltaR.h"

#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>

#include <FWCore/ParameterSet/interface/FileInPath.h>

#define DEBUG false

JetMETExtractor::JetMETExtractor(const std::string& name, const edm::ParameterSet& config)
: BaseExtractor(name, config)
{

  m_jetsTag = config.getParameter<edm::InputTag>("input_jets");
  m_metTag = config.getParameter<edm::InputTag>("input_met"); 
  m_rawMetTag = config.getParameter<edm::InputTag>("input_raw_met");
  m_particleFlowTag = config.getParameter<edm::InputTag>("pf_candidates");
  m_primaryVerticesTag = config.getParameter<edm::InputTag>("vertices");
  m_rhoTag = config.getParameter<edm::InputTag>("rho");

  mCorrectJets = config.getUntrackedParameter<bool>("redoJetCorrection", false);
  mUseGlobalTagForJEC = config.getUntrackedParameter<bool>("useGlobalTagForJEC", true);
  mUseType1Fix = config.getUntrackedParameter<bool>("useType1Fix", false);
  mUseGlobalTagForType1Fix = config.getUntrackedParameter<bool>("useGlobalTagForType1Fix", true);

  if (!mUseGlobalTagForJEC) {
    mJecPayload =  config.getUntrackedParameter<std::string>("jecPayload");
    mJecJetAlgo =  config.getUntrackedParameter<std::string>("jecJetAlgo");
    if (mJecPayload.length() > 0) {
      mJecPayload = edm::FileInPath(mJecPayload).fullPath();
    } else {
      std::cout << "WARNING! No JecPayload file found. Use the global tag instead for JEC" << std::endl;
      mUseGlobalTagForJEC = true;
    }
  }

  if (mUseType1Fix) {
    if (!mUseGlobalTagForType1Fix) {
      mJecPayload_L1ForType1Fix =  config.getUntrackedParameter<std::string>("jecPayload_L1ForType1Fix");
      mJecJetAlgo =  config.getUntrackedParameter<std::string>("jecJetAlgo");
      if (mJecPayload_L1ForType1Fix.length() > 0) {
        mJecPayload_L1ForType1Fix = edm::FileInPath(mJecPayload_L1ForType1Fix).fullPath();
      } else {
        std::cout << "WARNING! No JecPayload_L1ForType1Fix files found. Use the nominal Type I for MET correction" << std::endl;
        mUseGlobalTagForType1Fix = true;
      }
    }
  }

  mTxtCorrector = nullptr;
  mTxtCorrector_L1ForType1Fix = nullptr;

  if (mCorrectJets)
    mJetCorrectorLabel = config.getParameter<std::string>("jetCorrectorLabel");
  if (mUseType1Fix)
    mJetCorrectorLabelForType1Fix = config.getParameter<std::string>("jetCorrectorLabelForType1Fix");
  mDoJER = config.getUntrackedParameter<bool>("doJER", true);
  mJERSign = 0;
  if (mDoJER) {
    mJERSign = config.getUntrackedParameter<int>("jerSign", 0);
  }
  mDoLooseJetID = config.getUntrackedParameter<bool>("doLooseJetID", true);

  mJESSign = config.getUntrackedParameter<int>("jesSign", 0);

  if (mJERSign != 0 && mJERSign != -1 && mJERSign != 1)
    throw edm::Exception(edm::errors::LogicError) << "jerSign must be 0 for nominal correction, -1 for 1-sigma down correction, or 1 for 1-sigma up correction";

  if (mJESSign != 0 && mJESSign != -1 && mJESSign != 1)
    throw edm::Exception(edm::errors::LogicError) << "jesSign must be 0 for nominal correction, -1 for 1-sigma down correction, or 1 for 1-sigma up correction";

  mCorrectSysShiftMet = config.getUntrackedParameter<bool>("redoMetPhiCorrection", false);
  mRedoTypeI  = config.getUntrackedParameter<bool>("redoMetTypeICorrection", false);
  mSaveUnclusteredParticles = config.getUntrackedParameter<bool>("saveUnclusteredParticles", false);

  mJecFilename = config.getUntrackedParameter<std::string>("jes_uncertainties_file", "");
  if (mJESSign != 0 && mJecFilename.length() > 0) {
    mJecFilename = edm::FileInPath(mJecFilename).fullPath();
  }

#if DEBUG
  std::cout << "##########" << std::endl;
  std::cout << "JetMET extractor summary" << std::endl;
  std::cout << "Jet name: " << name << "; MET name: " << met_name << std::endl;
  std::cout << "Redo JEC: " << mCorrectJets << ((mCorrectJets) ? ("; " + mJetCorrectorLabel) : "") << std::endl;
  std::cout << "Do JER: " << mDoJER << "; JER sign: " << mJERSign << std::endl;
#endif

  // Set everything to 0
  m_jet_lorentzvector = new TClonesArray("TLorentzVector");
  m_genjet_lorentzvector = new TClonesArray("TLorentzVector");
  m_rawjet_lorentzvector = new TClonesArray("TLorentzVector");
  m_met_lorentzvector = new TClonesArray("TLorentzVector");
  m_unclustered_particle_lorentzvector = new TClonesArray("TLorentzVector");
  m_genmet_lorentzvector = new TClonesArray("TLorentzVector");
  m_scaleFactors.setWriteMode();

  reset();

  std::string jetsTreeName = config.getParameter<std::string>("tree_name_jets");
  std::string metTreeName = config.getParameter<std::string>("tree_name_met");

  // Tree definition

  m_tree_jet = NULL;
  m_tree_jet     = new TTree(jetsTreeName.c_str(), "PAT PF jet info");  
  m_tree_jet->SetAutoSave(0);
  m_tree_jet->Branch("n_jets",  &m_size,   "n_jets/i");  
  m_tree_jet->Branch("rho", &m_rho, "rho/F");
  m_tree_jet->Branch("jet_4vector","TClonesArray",&m_jet_lorentzvector, 5000, 0);
  m_tree_jet->Branch("genjet_4vector","TClonesArray",&m_genjet_lorentzvector, 5000, 0);
  m_tree_jet->Branch("rawjet_4vector","TClonesArray",&m_rawjet_lorentzvector, 5000, 0);
  m_tree_jet->Branch("jet_vx",  &m_jet_vx,   "jet_vx[n_jets]/F");  
  m_tree_jet->Branch("jet_vy",  &m_jet_vy,   "jet_vy[n_jets]/F");  
  m_tree_jet->Branch("jet_vz",  &m_jet_vz,   "jet_vz[n_jets]/F"); 
  m_tree_jet->Branch("jet_chmult",        &m_jet_chmult,       "jet_chmult[n_jets]/I");
  m_tree_jet->Branch("jet_chmuEfrac",     &m_jet_chmuEfrac,    "jet_chmuEfrac[n_jets]/F");
  m_tree_jet->Branch("jet_chemEfrac",     &m_jet_chemEfrac,    "jet_chemEfrac[n_jets]/F");
  m_tree_jet->Branch("jet_chhadEfrac",    &m_jet_chhadEfrac,   "jet_chhadEfrac[n_jets]/F");
  m_tree_jet->Branch("jet_nemEfrac",      &m_jet_nemEfrac,     "jet_nemEfrac[n_jets]/F");
  m_tree_jet->Branch("jet_nhadEfrac",     &m_jet_nhadEfrac,    "jet_nhadEfrac[n_jets]/F");
  m_tree_jet->Branch("jet_isPFJetLoose",  &m_jet_isPFJetLoose, "jet_isPFJetLoose[n_jets]/I");


  m_tree_jet->Branch("jet_puJetFullDiscriminant", &m_jet_puJetFullDiscriminant, "jet_puJetFullDiscriminant[n_jets]/F");
  m_tree_jet->Branch("jet_puJetFullId", &m_jet_puJetFullId, "jet_puJetFullId[n_jets]/I");
  m_tree_jet->Branch("jet_puJetCutBasedId", &m_jet_puJetCutBasedId, "jet_puJetCutBasedId[n_jets]/I");

  // 2012 b-tag algo
  m_tree_jet->Branch("jet_btag_jetProb", &m_jet_btag_jetProb, "jet_btag_jetProb[n_jets]/F");
  m_tree_jet->Branch("jet_btag_TCHP", &m_jet_btag_TCHP, "jet_btag_TCHP[n_jets]/F");
  m_tree_jet->Branch("jet_btag_CSV", &m_jet_btag_CSV, "jet_btag_CSV[n_jets]/F");
  // 2015
  m_tree_jet->Branch("jet_btag_CSVInclusive", &m_jet_btag_CSVInclusive, "jet_btag_CSVInclusive[n_jets]/F");

  m_tree_jet->Branch("jet_qgtag_likelihood", &m_jet_qgtag_likelihood, "jet_qgtag_likelihood[n_jets]/F");

  m_tree_jet->Branch("jet_mcParticleIndex",&m_jet_MCIndex,"jet_mcParticleIndex[n_jets]/I");

  m_tree_jet->Branch("jet_algo_parton_flavor", &m_jet_algo_parton_flavor, "jet_algo_parton_flavor[n_jets]/I");
  m_tree_jet->Branch("jet_physics_parton_pdgid", &m_jet_physics_parton_pdgid, "jet_physics_parton_pdgid[n_jets]/I");

  m_tree_jet->Branch("jet_scaleFactor", &m_scaleFactors.getBackingArray());
  m_tree_jet->Branch("jet_uncertainty_correctionFactor", &m_jet_uncertainty_correctionFactor, "jet_uncertainty_correctionFactor[n_jets]/F");

  m_tree_met = NULL;
  m_tree_met = new TTree(metTreeName.c_str(), "PAT PF MET info");  
  m_tree_met->SetAutoSave(0);
  m_tree_met->Branch("met_4vector","TClonesArray",&m_met_lorentzvector, 1000, 0);
  m_tree_met->Branch("sumEt", &m_met_sumEt, "sumEt/F");  
  m_tree_met->Branch("unclustered_particle_4vector","TClonesArray",&m_unclustered_particle_lorentzvector, 1000, 0);
  m_tree_met->Branch("genmet_4vector","TClonesArray",&m_genmet_lorentzvector, 1000, 0);
  
  
  // Mark that the extractor has been constructed properly
  setHealthy(true);
}

void JetMETExtractor::doConsumes(edm::ConsumesCollector&& collector) {
  BaseExtractor::doConsumes(std::forward<edm::ConsumesCollector>(collector));

  m_jetsToken = collector.consumes<pat::JetCollection>(m_jetsTag);
  m_metToken = collector.consumes<pat::METCollection>(m_metTag);

  if (mCorrectJets || mRedoTypeI) {
    m_rawMetToken = collector.consumes<pat::METCollection>(m_rawMetTag);
  }

  if (mSaveUnclusteredParticles) {
    m_particleFlowToken = collector.consumes<reco::PFCandidateCollection>(m_particleFlowTag);
  }

  m_primaryVerticesToken = collector.consumes<reco::VertexCollection>(m_primaryVerticesTag);
  m_rhoToken = collector.consumes<double>(m_rhoTag);
}

void JetMETExtractor::beginJob(bool isInAnalysisMode) {

  if (isInAnalysisMode)
    return;

  if (!mUseGlobalTagForJEC) {
    mTxtCorrector = makeFactorizedJetCorrectorFromXML(mJecPayload, mJecJetAlgo, m_isMC);
    std::cout << "Using text files for JEC" << std::endl;
  } else {
    std::cout << "Using global tag for JEC" << std::endl;
  }
  if (mUseType1Fix) {
    std::cout << "Do type 1 fix for MET" << std::endl;
    if (!mUseGlobalTagForType1Fix) {
      mTxtCorrector_L1ForType1Fix = makeFactorizedJetCorrectorFromXML(mJecPayload_L1ForType1Fix, mJecJetAlgo, m_isMC);
      std::cout << "Using text files for type 1 fix" << std::endl;
    } else {
      std::cout << "Using global for type 1 fix" << std::endl;
    }
  }
}


JetMETExtractor::JetMETExtractor(const std::string& name, const edm::ParameterSet& config, TFile *a_file)
: BaseExtractor(name, config, a_file)
{
  std::cout << "JetMETExtractor objet is retrieved" << std::endl;
  m_file = a_file;

  // Tree definition
  std::string jetsTreeName = config.getParameter<std::string>("tree_name_jets");
  std::string metTreeName = config.getParameter<std::string>("tree_name_met");

  m_tree_jet = dynamic_cast<TTree*>(a_file->Get(jetsTreeName.c_str()));


  if (m_tree_jet) {
    m_jet_lorentzvector = new TClonesArray("TLorentzVector");
    m_genjet_lorentzvector = new TClonesArray("TLorentzVector");
    m_rawjet_lorentzvector = new TClonesArray("TLorentzVector");

    if (m_tree_jet->FindBranch("n_jets")) 
      m_tree_jet->SetBranchAddress("n_jets",            &m_size);
    if (m_tree_jet->FindBranch("rho"))
      m_tree_jet->SetBranchAddress("rho",            &m_rho);
    if (m_tree_jet->FindBranch("jet_4vector")) 
      m_tree_jet->SetBranchAddress("jet_4vector",       &m_jet_lorentzvector);
    if (m_tree_jet->FindBranch("genjet_4vector")) 
      m_tree_jet->SetBranchAddress("genjet_4vector",       &m_genjet_lorentzvector);
    if (m_tree_jet->FindBranch("rawjet_4vector")) 
      m_tree_jet->SetBranchAddress("rawjet_4vector",       &m_rawjet_lorentzvector);
    if (m_tree_jet->FindBranch("jet_vx")) 
      m_tree_jet->SetBranchAddress("jet_vx",            &m_jet_vx);
    if (m_tree_jet->FindBranch("jet_vy")) 
      m_tree_jet->SetBranchAddress("jet_vy",            &m_jet_vy);
    if (m_tree_jet->FindBranch("jet_vz")) 
      m_tree_jet->SetBranchAddress("jet_vz",            &m_jet_vz);
    if (m_tree_jet->FindBranch("jet_uncertainty_correctionFactor"))
      m_tree_jet->SetBranchAddress("jet_uncertainty_correctionFactor",            &m_jet_uncertainty_correctionFactor);
    if (m_tree_jet->FindBranch("jet_chmult")) 
      m_tree_jet->SetBranchAddress("jet_chmult",        &m_jet_chmult);
    if (m_tree_jet->FindBranch("jet_chmuEfrac")) 
      m_tree_jet->SetBranchAddress("jet_chmuEfrac",     &m_jet_chmuEfrac);
    if (m_tree_jet->FindBranch("jet_chemEfrac")) 
      m_tree_jet->SetBranchAddress("jet_chemEfrac",     &m_jet_chemEfrac);
    if (m_tree_jet->FindBranch("jet_chhadEfrac")) 
      m_tree_jet->SetBranchAddress("jet_chhadEfrac",    &m_jet_chhadEfrac);
    if (m_tree_jet->FindBranch("jet_nemEfrac")) 
      m_tree_jet->SetBranchAddress("jet_nemEfrac",      &m_jet_nemEfrac);
    if (m_tree_jet->FindBranch("jet_nhadEfrac")) 
      m_tree_jet->SetBranchAddress("jet_nhadEfrac",     &m_jet_nhadEfrac);
    if (m_tree_jet->FindBranch("jet_isPFJetLoose")) 
      m_tree_jet->SetBranchAddress("jet_isPFJetLoose",  &m_jet_isPFJetLoose);

    if (m_tree_jet->FindBranch("jet_puJetFullDiscriminant"))
      m_tree_jet->SetBranchAddress("jet_puJetFullDiscriminant",  &m_jet_puJetFullDiscriminant);
    if (m_tree_jet->FindBranch("jet_puJetFullId"))
      m_tree_jet->SetBranchAddress("jet_puJetFullId",  &m_jet_puJetFullId);
    if (m_tree_jet->FindBranch("jet_puJetCutBasedId"))
      m_tree_jet->SetBranchAddress("jet_puJetCutBasedId",  &m_jet_puJetCutBasedId);

    if (m_tree_jet->FindBranch("jet_btag_jetProb"))
      m_tree_jet->SetBranchAddress("jet_btag_jetProb",  &m_jet_btag_jetProb);
    if (m_tree_jet->FindBranch("jet_btag_TCHP"))
      m_tree_jet->SetBranchAddress("jet_btag_TCHP",     &m_jet_btag_TCHP);
    if (m_tree_jet->FindBranch("jet_btag_CSV")) 
      m_tree_jet->SetBranchAddress("jet_btag_CSV",      &m_jet_btag_CSV);
    if (m_tree_jet->FindBranch("jet_btag_CSVInclusive")) 
      m_tree_jet->SetBranchAddress("jet_btag_CSVInclusive", &m_jet_btag_CSVInclusive);

    if (m_tree_jet->FindBranch("jet_qgtag_likelihood"))
      m_tree_jet->SetBranchAddress("jet_qgtag_likelihood",      &m_jet_qgtag_likelihood);

    if (m_tree_jet->FindBranch("jet_algo_parton_flavor"))
      m_tree_jet->SetBranchAddress("jet_algo_parton_flavor", &m_jet_algo_parton_flavor);

    if (m_tree_jet->FindBranch("jet_physics_parton_pdgid"))
      m_tree_jet->SetBranchAddress("jet_physics_parton_pdgid", &m_jet_physics_parton_pdgid);

    if (m_tree_jet->FindBranch("jet_mcParticleIndex")) 
      m_tree_jet->SetBranchAddress("jet_mcParticleIndex",&m_jet_MCIndex);

    if (m_tree_jet->FindBranch("jet_scaleFactor"))
      m_tree_jet->SetBranchAddress("jet_scaleFactor", &m_scaleFactors.getBackingArray());
  }

  m_tree_met = dynamic_cast<TTree*>(a_file->Get(metTreeName.c_str()));

  if (!m_tree_met)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_met_lorentzvector = new TClonesArray("TLorentzVector");
  m_genmet_lorentzvector = new TClonesArray("TLorentzVector");

  if (m_tree_met->FindBranch("met_4vector"))
    m_tree_met->SetBranchAddress("met_4vector", &m_met_lorentzvector);

  if (m_tree_met->FindBranch("sumEt")) 
      m_tree_met->SetBranchAddress("sumEt", &m_met_sumEt);

  m_unclustered_particle_lorentzvector = new TClonesArray("TLorentzVector");

  if (m_tree_met->FindBranch("unclustered_particle_4vector"))
    m_tree_met->SetBranchAddress("unclustered_particle_4vector", &m_unclustered_particle_lorentzvector);

  if (m_tree_met->FindBranch("genmet_4vector"))
    m_tree_met->SetBranchAddress("genmet_4vector",&m_genmet_lorentzvector);
  
  
  // Mark that the extractor has been constructed properly
  setHealthy(true);
}


JetMETExtractor::~JetMETExtractor()
{
  delete mTxtCorrector;
  delete mTxtCorrector_L1ForType1Fix;
}


bool JetMETExtractor::isPFJetLoose(const pat::Jet& jet)
{
  // Taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetID

  if (! jet.isPFJet())
    return false;

  // Jet ID works only with uncorrected jet. *EnergyFraction functions take care of that all alone
  bool isValid = jet.neutralHadronEnergyFraction() < 0.99;
  isValid &= jet.neutralEmEnergyFraction() < 0.99;
  isValid &= jet.numberOfDaughters() > 1;
  if (fabs(jet.eta()) < 2.4) {
    isValid &= jet.chargedHadronEnergyFraction() > 0.;
    isValid &= jet.chargedMultiplicity() > 0;
    isValid &= jet.chargedEmEnergyFraction() < 0.99;
  }

  return isValid;
}

//
// Method filling the main particle tree
//

void JetMETExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* m_MC) 
{
  edm::Handle<pat::JetCollection>  jetHandle;
  event.getByToken(m_jetsToken, jetHandle);
  pat::JetCollection p_jets = *jetHandle;

  reset();
  m_size = 0;

  if (mJESSign != 0) {
    if (! jecUncertainty.get()) {
      if (mJecFilename.length() > 0) {
        std::cout << "Reading JES uncertainties from '" << mJecFilename << "'" << std::endl;
        jecUncertainty.reset(new JetCorrectionUncertainty(mJecFilename));
      } else {
        std::cout << "Reading JES uncertainties from Global Tag" << std::endl;
        edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
        iSetup.get<JetCorrectionsRecord>().get("AK5PFchs",JetCorParColl);
        JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
        jecUncertainty.reset(new JetCorrectionUncertainty(JetCorPar));
      }
    }
  }

  if (mCorrectJets) {
    correctJets(p_jets, event, iSetup);
  } else {
    extractRawJets(p_jets);
  }

  edm::Handle<pat::METCollection> metHandle;
  event.getByToken(m_metToken, metHandle);
  pat::MET MET = (*metHandle).at(0);

  // If we are redoing jet correction, or if the user forces it, recompute TypeI MEt corrections
  if (mCorrectJets || mRedoTypeI) {
    // Raw MET
    edm::Handle<pat::METCollection> rawMetHandle;
    event.getByToken(m_rawMetToken, rawMetHandle);

    if (rawMetHandle.isValid()) {
      const pat::MET& rawMet = rawMetHandle->back();
      correctMETWithTypeI(rawMet, MET, p_jets, event, iSetup);
    }
  }

  //do Jets MET resolution corrections
  if (m_MC && mDoJER)
    correctJetsMETresolution(p_jets, MET);

  //do MET SysShift corrections
  if (mCorrectSysShiftMet)
    correctMETWithSysShift(event, MET);

  // JES systematics
  if (mJESSign != 0)
    doJESSystematics(p_jets, MET);


  for (unsigned int i = 0; i < p_jets.size(); ++i)
  {
    const pat::Jet& rawJet = *(p_jets.at(i).userData<pat::Jet>("rawJet"));
#if DEBUG
      if (! mCorrectJets) {
        // Test if jet id works the same on our stored raw jets and the jet itself
        if (isPFJetLoose(rawJet) != isPFJetLoose(p_jets.at(i))) {
          std::cout << "Error: there's something wrong with your jets!" << std::endl;
        }
      }
#endif
    if(mDoLooseJetID) {
      if (! isPFJetLoose(rawJet)) 
        continue;
    }

      
    pat::JetRef jetRef(jetHandle, i);    
    JetMETExtractor::writeInfo(event, iSetup, p_jets.at(i), m_size, jetRef); 


    // Retrieve QGTag info
    float qgLikelihood = p_jets.at(i).userFloat("QGTaggerPFlow:qgLikelihood");
    m_jet_qgtag_likelihood[i] = qgLikelihood;


    if (m_MC)
      doMCMatch(p_jets.at(i), event, m_MC, m_size);

    m_size++;
  }

  writeInfo(event, iSetup, MET, 0);

  if (mSaveUnclusteredParticles) {
    edm::Handle<reco::PFCandidateCollection>  PFParticleHandle;
    event.getByToken(m_particleFlowToken, PFParticleHandle);
    reco::PFCandidateCollection p_PFParticle = *PFParticleHandle;

    int PFParticle_size = 0;

    for (unsigned int i = 0; i < p_PFParticle.size(); ++i)
    {
      reco::PFCandidate& mainPFPart = p_PFParticle[i];
      int pdgId = mainPFPart.pdgId();
      if (fabs(pdgId) == 11 || fabs(pdgId) == 13)
        continue;
      bool found = false;
      for (unsigned int j = 0; j < p_jets.size(); ++j) {
        pat::Jet& jet = p_jets[j];
        double mainDeltaR = deltaR(mainPFPart, jet);
        if (mainDeltaR < 0.5) {
          found = true;
          break;
        }
        if (found)
          break;
      }
      if (!found) {
        writeInfo(event, iSetup, mainPFPart, PFParticle_size);
        PFParticle_size ++;
      }
    }
  }


  fillTree();
}

void JetMETExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Jet& part, int index, const pat::JetRef& ref) 
{
  if (index>=m_jets_MAX) return;

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Writing jet #" << index << std::endl;
  std::cout << "Pt: " << part.pt() << "; Px / Py / Pz / E : " << part.px() << " / " << part.py() << " / " << part.pz() << " / " << part.energy() << std::endl;
  std::cout << "Eta: " << part.eta() << "; Phi : " << part.phi() << std::endl;
#endif

  new((*m_jet_lorentzvector)[index]) TLorentzVector(part.px(),part.py(),part.pz(),part.energy());
  
  if(part.genJet()) {
    new((*m_genjet_lorentzvector)[index]) TLorentzVector((part.genJet())->px(),(part.genJet())->py(),(part.genJet())->pz(),(part.genJet())->energy());
  }
  else {
    new((*m_genjet_lorentzvector)[index]) TLorentzVector(0.,0.,0.,0.);
  }

  const pat::Jet* rawJet = part.userData<pat::Jet>("rawJet");
  new((*m_rawjet_lorentzvector)[index]) TLorentzVector(rawJet->px(),rawJet->py(),rawJet->pz(),rawJet->energy());

  m_jet_vx[index]   = part.vx();
  m_jet_vy[index]   = part.vy();
  m_jet_vz[index]   = part.vz();

  if (part.isPFJet())
  {
    m_jet_chmult[index]        = part.chargedMultiplicity();
    m_jet_chmuEfrac[index]     = part.chargedMuEnergyFraction();
    m_jet_chemEfrac[index]     = part.chargedEmEnergyFraction();
    m_jet_chhadEfrac[index]    = part.chargedHadronEnergyFraction();
    m_jet_nemEfrac[index]      = part.neutralEmEnergyFraction();
    m_jet_nhadEfrac[index]     = part.neutralHadronEnergyFraction();
    
    m_jet_isPFJetLoose[index]  = int(isPFJetLoose(part));

    // PU Jet ID from JetToolBox
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetToolbox
    // for pattuples
    //m_jet_puJetFullDiscriminant[index] = part.hasUserFloat("pileupJetIdEvaluatorPFlow:fullDiscriminant") ? part.userFloat("pileupJetIdEvaluatorPFlow:fullDiscriminant") : -2;
    // for miniAOD
    m_jet_puJetFullDiscriminant[index] = part.hasUserFloat("pileupJetId:fullDiscriminant") ? part.userFloat("pileupJetId:fullDiscriminant") : -2;
    m_jet_puJetFullId[index]           = part.hasUserInt("pileupJetIdEvaluatorPFlow:fullId") ? part.userFloat("pileupJetIdEvaluatorPFlow:fullId") : -1;
    m_jet_puJetCutBasedId[index]       = part.hasUserInt("pileupJetIdEvaluatorPFlow:cutbasedId") ? part.userFloat("pileupJetIdEvaluatorPFlow:cutbasedId") : -1;

    // 2012: See https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
    m_jet_btag_jetProb[index]      = part.bDiscriminator("jetProbabilityBJetTags");
    m_jet_btag_TCHP[index]         = part.bDiscriminator("trackCountingHighPurBJetTags");
    m_jet_btag_CSV[index]          = part.bDiscriminator("combinedSecondaryVertexBJetTags");

    // 2015
    m_jet_btag_CSVInclusive[index] = part.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags");

    if (part.hasUserFloat("uncertainty_correctionFactor"))
      m_jet_uncertainty_correctionFactor[index] = part.userFloat("uncertainty_correctionFactor");
    else
      m_jet_uncertainty_correctionFactor[index] = 1.;
  }

  m_jet_algo_parton_flavor[index] = part.partonFlavour();
  m_jet_physics_parton_pdgid[index] = (part.genParton()) ? part.genParton()->pdgId() : 0;

  if (m_isMC) {
    int mcFlavor = abs(m_jet_algo_parton_flavor[index]);
    ScaleFactorService::Flavor flavor = ScaleFactorService::B;
    if (mcFlavor == 4) {
      flavor = ScaleFactorService::C;
    } else if ((mcFlavor <= 3) || (mcFlavor == 21)) {
      // If mcFlavor == 0, assume it's a light jet
      flavor = ScaleFactorService::LIGHT;
    }

    m_scaleFactors.push_back(ScaleFactorService::getInstance().getBTaggingScaleFactor(flavor, part.et(), part.eta()));
  }
}

void JetMETExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::MET& part, int index) 
{
  if (index > 1)
    return;

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Writing MET #" << index << std::endl;
  std::cout << "Pt: " << part.pt() << "; Px / Py / Pz : " << part.px() << " / " << part.py() << " / " << part.pz() << std::endl;
  std::cout << "Eta: " << part.eta() << "; Phi : " << part.phi() << std::endl;
#endif

  new((*m_met_lorentzvector)[index]) TLorentzVector(part.px(),part.py(),part.pz(),part.energy());

  if (part.genMET()) {
    const reco::GenMET* genmet = part.genMET();
    new((*m_genmet_lorentzvector)[index]) TLorentzVector(genmet->px(),genmet->py(),genmet->pz(),genmet->energy());  
  } else {
    new((*m_genmet_lorentzvector)[index]) TLorentzVector(0., 0., 0., 0.);
  }

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Writing GenMET #" << index << std::endl;
  std::cout << "Pt: " << genmet->pt() << "; Px / Pz / Pz : " << genmet->px() << " / " << genmet->py() << " / " << genmet->pz() << std::endl;
  std::cout << "Eta: " << genmet->eta() << "; Phi : " << genmet->phi() << std::endl;
#endif

}

void JetMETExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const reco::PFCandidate& part, int index)
{
#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Writing PF particle #" << index << std::endl;
  std::cout << "Pt: " << part.pt() << "; Px / Py / Pz / E : " << part.px() << " / " << part.py() << " / " << part.pz() << " / " << part.energy() << std::endl;
  std::cout << "Eta: " << part.eta() << "; Phi : " << part.phi() << std::endl;
#endif

  new((*m_unclustered_particle_lorentzvector)[index]) TLorentzVector(part.px(),part.py(),part.pz(),part.energy());
}

//
// Method getting the info from an input file
//

void JetMETExtractor::getInfo(int ievt) 
{
  if (m_tree_jet)
    m_tree_jet->GetEntry(ievt);

  if (m_tree_met)
    m_tree_met->GetEntry(ievt);
}


// Method initializing everything (to do for each event)

void JetMETExtractor::reset()
{
  m_size = 0;
  m_rho = 10000.;

  m_scaleFactors.clear();

  for (int i=0;i<m_jets_MAX;++i) 
  {
    m_jet_vx[i] = 0.;
    m_jet_vy[i] = 0.;
    m_jet_vz[i] = 0.;

    m_jet_uncertainty_correctionFactor[i] = 0.;

    m_jet_chmult[i] = 0;
    m_jet_chmuEfrac[i] = 0.;
    m_jet_chemEfrac[i] = 0.;
    m_jet_chhadEfrac[i] = 0.;
    m_jet_nemEfrac[i] = 0.;
    m_jet_nhadEfrac[i] = 0.;
    m_jet_isPFJetLoose[i] = 0.;
    m_jet_puJetFullDiscriminant[i] = -3.;
    m_jet_puJetFullId[i] = 0.;
    m_jet_puJetCutBasedId[i] = 0.;
    
    m_jet_btag_jetProb[i] = 0.;
    m_jet_btag_TCHP[i] = 0.;
    m_jet_btag_CSV[i] = 0.;
    m_jet_btag_CSVInclusive[i] = 0.;

    m_jet_MCIndex[i]    = -1;

    m_jet_algo_parton_flavor[i] = 0;
    m_jet_physics_parton_pdgid[i] = 0;
  }
  if (m_jet_lorentzvector)
    m_jet_lorentzvector->Clear();
  if (m_genjet_lorentzvector)
    m_genjet_lorentzvector->Clear();
  if (m_rawjet_lorentzvector)
    m_rawjet_lorentzvector->Clear();
  if (m_met_lorentzvector)
    m_met_lorentzvector->Clear();

  m_met_sumEt = 0.;

  if (m_unclustered_particle_lorentzvector)
    m_unclustered_particle_lorentzvector->Clear();

  if (m_genmet_lorentzvector)
    m_genmet_lorentzvector->Clear();
}


void JetMETExtractor::fillTree()
{
  if (m_tree_jet)
    m_tree_jet->Fill(); 

  if (m_tree_met)
    m_tree_met->Fill();
}

void JetMETExtractor::correctJets(pat::JetCollection& jets, const edm::Event& iEvent, const edm::EventSetup& iSetup) {

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Recompute jet energy corrections..." << std::endl;
#endif

  // Get Jet corrector
  const JetCorrector* globalTagCorrector = nullptr;

  if (mUseGlobalTagForJEC) {
    globalTagCorrector = JetCorrector::getJetCorrector(mJetCorrectorLabel, iSetup);
  }

  edm::Handle<reco::VertexCollection>  vertexHandle;
  iEvent.getByToken(m_primaryVerticesToken, vertexHandle);
  reco::VertexCollection vertices = *vertexHandle;

  edm::Handle<double> rhos;
  iEvent.getByToken(m_rhoToken, rhos);
  double rho = *rhos;

  // Correct jets
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    pat::Jet& jet = *it;

    pat::Jet rawJet = jet.correctedJet("Uncorrected");
    jet.addUserData("rawJet", rawJet, true); // Store raw jet inside the jet
    const pat::Jet L1Jet  = jet.correctedJet("L1FastJet");
    jet.addUserData("L1Jet", L1Jet, true); // Embed L1 corrected jet for TypeI correction

#if DEBUG
    std::cout << "---" << std::endl;
    std::cout << "Pt: " << jet.pt() << std::endl;
#endif

    double toRaw = jet.jecFactor("Uncorrected");
    jet.setP4(jet.p4() * toRaw); // jet is now a raw jet
#if DEBUG
    std::cout << "True raw pt: " << rawJet.pt() << std::endl;
    std::cout << "Raw pt: " << jet.pt() << std::endl;
#endif

    double corrections = 0.;
    if (mUseGlobalTagForJEC) {
      corrections = globalTagCorrector->correction(jet, iEvent, iSetup);
    } else {
      mTxtCorrector->setJetEta(jet.eta());
      mTxtCorrector->setJetPt(jet.pt());
      mTxtCorrector->setRho(rho);
      mTxtCorrector->setJetA(jet.jetArea());
      mTxtCorrector->setNPV(vertices.size());
      corrections = mTxtCorrector->getCorrection();
    }
       
    jet.scaleEnergy(corrections);
#if DEBUG
    std::cout << "Corrected pt: " << jet.pt() << std::endl;
#endif
  }

  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorter);
}

//from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution

double JetMETExtractor::getResCorrFactor(const pat::Jet& jet) {

  const auto& getFactor = [this] (float nominal, float down, float up) -> float {
    switch (mJERSign) {
      case 0:
        return nominal;

      case 1:
        return up;

      case -1:
        return down;

      default:
        throw new edm::Exception(edm::errors::LogicError, "Invalid JER sign: " + mJERSign);
    }
  };

  // JER Scaling factors and Uncertainty for 13 TeV (2015)
  if (fabs(jet.eta()) > 0. && fabs(jet.eta()) <= 0.8) {
    return getFactor(1.061, 1.038, 1.084);
  } else if (fabs(jet.eta()) > 0.8 && fabs(jet.eta()) <= 1.3) {
    return getFactor(1.088, 1.059, 1.117);
  } else if (fabs(jet.eta()) > 1.3 && fabs(jet.eta()) <= 1.9) {
    return getFactor(1.106, 1.076, 1.136);
  } else if (fabs(jet.eta()) > 1.9 && fabs(jet.eta()) <= 2.5) {
    return getFactor(1.126, 1.032, 1.220);
  } else if (fabs(jet.eta()) > 2.5 && fabs(jet.eta()) <= 3.0) {
    return getFactor(1.343, 1.220, 1.466);
  } else if (fabs(jet.eta()) > 3.0 && fabs(jet.eta()) <= 3.2) {
    return getFactor(1.303, 1.192, 1.414);
  } else if (fabs(jet.eta()) > 3.2 && fabs(jet.eta()) <= 5.0) {
    return getFactor(1.320, 1.034, 1.606);
  }


  // JER Scaling factors and Uncertainty for 8 TeV (2012)
/*  if (fabs(jet.eta()) > 0. && fabs(jet.eta()) <= 0.5) {*/
    //return getFactor(1.079, 1.053, 1.105);
  //} else if (fabs(jet.eta()) > 0.5 && fabs(jet.eta()) <= 1.1) {
    //return getFactor(1.099, 1.071, 1.127);
  //} else if (fabs(jet.eta()) > 1.1 && fabs(jet.eta()) <= 1.7) {
    //return getFactor(1.121, 1.092, 1.150);
  //} else if (fabs(jet.eta()) > 1.7 && fabs(jet.eta()) <= 2.3) {
    //return getFactor(1.208, 1.162, 1.254);
  //} else if (fabs(jet.eta()) > 2.3 && fabs(jet.eta()) <= 2.8) {
    //return getFactor(1.254, 1.192, 1.316);
  //} else if (fabs(jet.eta()) > 2.8 && fabs(jet.eta()) <= 3.2) {
    //return getFactor(1.395, 1.332, 1.458);
  //} else if (fabs(jet.eta()) > 3.2 && fabs(jet.eta()) <= 5.0) {
    //return getFactor(1.056, 0.865, 1.247);
  /*}*/

  return 1.;
}

void JetMETExtractor::correctJetsMETresolution(pat::JetCollection& jets, pat::MET& met) {

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Doing jet resolution smearing" << std::endl;
#endif
  
  double correctedMetPx = met.px(); 
  double correctedMetPy = met.py(); 
  // Correct jets
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    pat::Jet& jet = *it;

    if (jet.pt() > 10) {

      const pat::Jet* rawJet = jet.userData<pat::Jet>("rawJet");

#if DEBUG
      std::cout << "---" << std::endl;
      std::cout << "Pt: " << jet.pt() << std::endl;
#endif

      // resolution corection as in https://twiki.cern.ch/twiki/bin/view/CMS/TWikiTopRefSyst#Jet_energy_resolution
      
      if (jet.genJet() == NULL)
        continue;

      double genjet_pt = jet.genJet()->pt();
      double jet_pt = jet.pt();
      double rescorr = getResCorrFactor(jet) - 1;
      double deltapt = (jet_pt - genjet_pt) * rescorr; 
      double scalefac = (jet_pt + deltapt) / jet_pt;
      if (scalefac <= 0)
        continue;
      
      correctedMetPx += (rawJet->px());
      correctedMetPy += (rawJet->py());
      m_met_sumEt -= (rawJet->pt());

      jet.scaleEnergy(scalefac);
      //rawJet->scaleEnergy(ptscale); 

#if DEBUG
      std::cout << "Corrected pt: " << jet.pt() << std::endl;
#endif
    
      correctedMetPx -= (rawJet->px()*scalefac);
      correctedMetPy -= (rawJet->py()*scalefac);
      m_met_sumEt += (rawJet->pt()*scalefac);
    }
  }

  
  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);
    
  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));
    

  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorter);
}
    

void JetMETExtractor::correctMETWithTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets, const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  double deltaPx = 0., deltaPy = 0., deltaPt = 0.;
#if DEBUG
    std::cout << "---" << std::endl;
    std::cout << "Computing TypeI correction" << std::endl;
    std::cout << "MET Raw et: " << rawMet.et() << std::endl;
    std::cout << "PAT corrected MET et: " << met.et() << std::endl;
#endif

  edm::Handle<reco::VertexCollection>  vertexHandle;
  iEvent.getByToken(m_primaryVerticesToken, vertexHandle);
  reco::VertexCollection vertices = *vertexHandle;

  edm::Handle<double> rhos;
  iEvent.getByToken(m_rhoToken, rhos);
  double rho = *rhos;
  m_rho = *rhos;

  // See https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=174324 slide 4
  // and http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h?revision=1.8&view=markup
  for (pat::JetCollection::const_iterator it = jets.begin(); it != jets.end(); ++it) {

    pat::Jet jet;
    pat::Jet rawJet = *(it->userData<pat::Jet>("rawJet"));
    pat::Jet L1Jet;
    if (! mUseType1Fix) {
      jet = *it;
      L1Jet  = *(it->userData<pat::Jet>("L1Jet"));
    } else {
      jet = rawJet;
      L1Jet = rawJet;

      const JetCorrector* globalTagCorrector = nullptr;
      const JetCorrector* globalTagCorrectorForType1Fix = nullptr;

      double jet_corrections = 0.;
      double L1Jet_corrections = 0.;
      if (mUseGlobalTagForType1Fix) {
        globalTagCorrector = JetCorrector::getJetCorrector(mJetCorrectorLabel, iSetup);
        globalTagCorrectorForType1Fix = JetCorrector::getJetCorrector(mJetCorrectorLabelForType1Fix, iSetup);
        jet_corrections = globalTagCorrector->correction(jet, iEvent, iSetup);
        jet_corrections = globalTagCorrectorForType1Fix->correction(L1Jet, iEvent, iSetup);
      } else {
        mTxtCorrector->setJetEta(jet.eta());
        mTxtCorrector->setJetPt(jet.pt());
        mTxtCorrector->setRho(rho);
        mTxtCorrector->setJetA(rawJet.jetArea());
        mTxtCorrector->setNPV(vertices.size());
        jet_corrections = mTxtCorrector->getCorrection();

        mTxtCorrector_L1ForType1Fix->setJetEta(L1Jet.eta());
        mTxtCorrector_L1ForType1Fix->setJetPt(L1Jet.pt());
        mTxtCorrector_L1ForType1Fix->setRho(rho);
        mTxtCorrector_L1ForType1Fix->setJetA(rawJet.jetArea());
        mTxtCorrector_L1ForType1Fix->setNPV(vertices.size());
        L1Jet_corrections = mTxtCorrector_L1ForType1Fix->getCorrection();
      }

      jet.scaleEnergy(jet_corrections);
      L1Jet.scaleEnergy(L1Jet_corrections);

    }

    if (jet.pt() > 10) {

      double emEnergyFraction = rawJet.chargedEmEnergyFraction() + rawJet.neutralEmEnergyFraction();
      if (emEnergyFraction > 0.90)
        continue;

      //reco::Candidate::LorentzVector rawJetP4 = rawJet->p4();
      reco::Candidate::LorentzVector L1JetP4  = L1Jet.p4();

      // Skip muons
      /*std::vector<reco::PFCandidatePtr> cands = rawJet->getPFConstituents();
        for (std::vector<reco::PFCandidatePtr>::const_iterator cand = cands.begin(); cand != cands.end(); ++cand) {
        if ((*cand)->muonRef().isNonnull() && skipMuonSelection(*(*cand)->muonRef())) {
        reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
        rawJetP4 -= muonP4;
        }
        }*/


      deltaPx += (jet.px() - L1JetP4.px());
      deltaPy += (jet.py() - L1JetP4.py());
      deltaPt += (jet.pt() - L1JetP4.pt());
    }
  }

  double correctedMetPx = rawMet.px() - deltaPx;
  double correctedMetPy = rawMet.py() - deltaPy;
  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);

  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));
  m_met_sumEt = rawMet.sumEt() + deltaPt;

#if DEBUG
    std::cout << "Handmade corrected MET et: " << met.et() << std::endl;
    std::cout << "Raw MET sumEt: " << rawMet.sumEt() << std::endl;
    std::cout << "Handmade corrected MET sumEt: " << m_met_sumEt() << std::endl;
#endif
}


//from JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi, check updates
//updated to MET systematic shift corrections to 2012 ABCD ReReco data + new Summer'13 JEC
double JetMETExtractor::getSysShifCorrFactorX(const int Nvtx){
  if (m_isMC) return -(+1.62861e-01 - 2.38517e-02*Nvtx);
  else        return -(+4.83642e-02 + 2.48870e-01*Nvtx);
}

double JetMETExtractor::getSysShifCorrFactorY(const int Nvtx){
  if (m_isMC) return -(+3.60860e-01 - 1.30335e-01*Nvtx);
  else        return -(-1.50135e-01 - 8.27917e-02*Nvtx);
}

void JetMETExtractor::correctMETWithSysShift(const edm::Event& event, pat::MET& met) {

  edm::Handle<reco::VertexCollection>  vertexHandle;
  event.getByToken(m_primaryVerticesToken, vertexHandle);
  reco::VertexCollection vertices = *vertexHandle;

  int Nvtx=0;
  for (reco::VertexCollection::iterator it = vertices.begin(); it != vertices.end(); ++it)  {
    reco::Vertex& vertex= *it;
    //cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
    if ( ! vertex.isValid() )  continue;
    if ( vertex.ndof() <  4 ) continue;
    if ( vertex.chi2() <= 0 ) continue;
    if ( vertex.tracksSize() <= 0 ) continue;
    if ( fabs(vertex.z())   >= 24 ) continue;
    if ( fabs(vertex.position().Rho()) >= 2 ) continue;

    Nvtx++;
  }


#if DEBUG
  std::cout << "Correcting MET phi shift" << std::endl;
  std::cout << "MET et without SysShift corrections: " << met.et() << std::endl;
#endif

  double correctedMetPx = met.px() + getSysShifCorrFactorX(Nvtx);
  double correctedMetPy = met.py() + getSysShifCorrFactorY(Nvtx);
  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);

  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));

#if DEBUG
    std::cout << "MET et with SysShift corrections: " << met.et() << std::endl;
#endif

}

void JetMETExtractor::doJESSystematics(pat::JetCollection& jets, pat::MET& met) {
#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "JES systematic" << std::endl;
#endif

  double correctedMetPx = met.px();
  double correctedMetPy = met.py();

  // Correct jets
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    pat::Jet& jet = *it;

    jecUncertainty->setJetEta(jet.eta());
    jecUncertainty->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt

    double uncertainty = (mJESSign == 1) ? fabs(jecUncertainty->getUncertainty(true)) : fabs(jecUncertainty->getUncertainty(false));
    double signedCorrection = mJESSign * uncertainty;

    double scaleFactor = 1. + signedCorrection;

    correctedMetPx += (jet.px());
    correctedMetPy += (jet.py());
    m_met_sumEt -= (jet.pt());

#if DEBUG
    std::cout << "---" << std::endl;
    std::cout << "Pt before JES uncertainty: " << jet.pt() << std::endl;
#endif

    jet.addUserFloat("uncertainty_correctionFactor", scaleFactor);
    jet.scaleEnergy(scaleFactor);

#if DEBUG
    std::cout << "Pt after JES uncertainty: " << jet.pt() << std::endl;
#endif

    correctedMetPx -= (jet.px());
    correctedMetPy -= (jet.py());
    m_met_sumEt += (jet.pt());
  }

  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "MET Et before JES uncertainty: " << met.et() << std::endl;
#endif

  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));

#if DEBUG
  std::cout << "MET Et after JES uncertainty: " << met.et() << std::endl;
#endif

  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorter);
}

void JetMETExtractor::extractRawJets(pat::JetCollection& jets) {

  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it) {
    pat::Jet& jet = *it;

    const pat::Jet rawJet = jet.correctedJet("Uncorrected");
    jet.addUserData("rawJet", rawJet, true);
    const pat::Jet L1Jet  = jet.correctedJet("L1FastJet");
    jet.addUserData("L1Jet", L1Jet, true); // Embed L1 corrected jet for TypeI correction
  }

}

void JetMETExtractor::setMETLorentzVector(int idx, float E, float Px, float Py, float Pz)
{
  new((*m_met_lorentzvector)[idx]) TLorentzVector(Px,Py,Pz,E);
}
