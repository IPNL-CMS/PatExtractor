#include "../interface/JetMETExtractor.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/VertexReco/interface/Vertex.h>

#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>

#include <FWCore/ParameterSet/interface/FileInPath.h>

#define DEBUG false

JetMETExtractor::JetMETExtractor(const std::string& name, const std::string& met_name, const edm::ParameterSet& config)
: BaseExtractor(name)
{
  if (! config.exists(name))
    throw edm::Exception(edm::errors::ConfigFileReadError) << "No edm::ParameterSet named " << name << " found";

  if (! config.exists(met_name))
    throw edm::Exception(edm::errors::ConfigFileReadError) << "No edm::ParameterSet named " << met_name << " found";

  const edm::ParameterSet& jetConfig = config.getParameter<edm::ParameterSet>(name);
  const edm::ParameterSet& metConfig = config.getParameter<edm::ParameterSet>(met_name);

  m_tag = jetConfig.getParameter<edm::InputTag>("input");
  m_metTag = metConfig.getParameter<edm::InputTag>("input"); 

  mCorrectJets = jetConfig.getUntrackedParameter<bool>("redoJetCorrection", false);
  if (mCorrectJets)
    mJetCorrectorLabel = jetConfig.getParameter<std::string>("jetCorrectorLabel");
  mDoJER = jetConfig.getUntrackedParameter<bool>("doJER", true);
  mJERSign = 0;
  if (mDoJER) {
    mJERSign = jetConfig.getUntrackedParameter<int>("jerSign", 0);
  }

  mJESSign = jetConfig.getUntrackedParameter<int>("jesSign", 0);

  if (mJERSign != 0 && mJERSign != -1 && mJERSign != 1)
    throw edm::Exception(edm::errors::LogicError) << "jerSign must be 0 for nominal correction, -1 for 1-sigma down correction, or 1 for 1-sigma up correction";

  if (mJESSign != 0 && mJESSign != -1 && mJESSign != 1)
    throw edm::Exception(edm::errors::LogicError) << "jesSign must be 0 for nominal correction, -1 for 1-sigma down correction, or 1 for 1-sigma up correction";

  mCorrectSysShiftMet = metConfig.getUntrackedParameter<bool>("redoMetPhiCorrection", false);
  mRedoTypeI  = metConfig.getUntrackedParameter<bool>("redoMetTypeICorrection", false);

  mJecFilename = jetConfig.getUntrackedParameter<std::string>("jes_uncertainties_file", "");
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
  m_met_lorentzvector = new TClonesArray("TLorentzVector");
  m_scaleFactors.setWriteMode();

  reset();

  // Tree definition

  m_tree_jet = NULL;
  m_tree_jet     = new TTree(name.c_str(), "PAT PF jet info");  
  m_tree_jet->Branch("n_jets",  &m_size,   "n_jets/i");  
  m_tree_jet->Branch("jet_4vector","TClonesArray",&m_jet_lorentzvector, 5000, 0);
  m_tree_jet->Branch("jet_vx",  &m_jet_vx,   "jet_vx[n_jets]/F");  
  m_tree_jet->Branch("jet_vy",  &m_jet_vy,   "jet_vy[n_jets]/F");  
  m_tree_jet->Branch("jet_vz",  &m_jet_vz,   "jet_vz[n_jets]/F"); 
  m_tree_jet->Branch("jet_chmult",        &m_jet_chmult,       "jet_chmult[n_jets]/I");
  m_tree_jet->Branch("jet_chmuEfrac",     &m_jet_chmuEfrac,    "jet_chmuEfrac[n_jets]/F");
  m_tree_jet->Branch("jet_chemEfrac",     &m_jet_chemEfrac,    "jet_chemEfrac[n_jets]/F");
  m_tree_jet->Branch("jet_chhadEfrac",    &m_jet_chhadEfrac,   "jet_chhadEfrac[n_jets]/F");
  m_tree_jet->Branch("jet_nemEfrac",      &m_jet_nemEfrac,     "jet_nemEfrac[n_jets]/F");
  m_tree_jet->Branch("jet_nhadEfrac",     &m_jet_nhadEfrac,    "jet_nhadEfrac[n_jets]/F");

  // 2012 b-tag algo
  m_tree_jet->Branch("jet_btag_jetProb",  &m_jet_btag_jetProb, "jet_btag_jetProb[n_jets]/F");
  m_tree_jet->Branch("jet_btag_TCHP",    &m_jet_btag_TCHP,   "jet_btag_TCHP[n_jets]/F");
  m_tree_jet->Branch("jet_btag_CSV",    &m_jet_btag_CSV,   "jet_btag_CSV[n_jets]/F");

  //m_tree_jet->Branch("jet_btag_BjetProb", &m_jet_btag_BjetProb,"jet_btag_BjetProb[n_jets]/F");
  //m_tree_jet->Branch("jet_btag_SSVHE",    &m_jet_btag_SSVHE,   "jet_btag_SSVHE[n_jets]/F");
  //m_tree_jet->Branch("jet_btag_SSVHP",    &m_jet_btag_SSVHP,   "jet_btag_SSVHP[n_jets]/F");
  //m_tree_jet->Branch("jet_btag_TCHE",    &m_jet_btag_TCHE,   "jet_btag_TCHE[n_jets]/F");
  m_tree_jet->Branch("jet_mcParticleIndex",&m_jet_MCIndex,"jet_mcParticleIndex[n_jets]/I");

  m_tree_jet->Branch("jet_scaleFactor", &m_scaleFactors.getBackingArray());

  m_tree_met = NULL;
  m_tree_met      = new TTree(met_name.c_str(), "PAT PF MET info");  
  m_tree_met->Branch("met_4vector","TClonesArray",&m_met_lorentzvector, 1000, 0);
}


JetMETExtractor::JetMETExtractor(const std::string& name, const std::string& met_name, TFile *a_file)
: BaseExtractor(name)
{
  std::cout << "JetMETExtractor objet is retrieved" << std::endl;
  m_file = a_file;

  // Tree definition
  m_OK = false;

  m_tree_jet = dynamic_cast<TTree*>(a_file->Get(name.c_str()));


  if (m_tree_jet) {
    m_jet_lorentzvector = new TClonesArray("TLorentzVector");

    if (m_tree_jet->FindBranch("n_jets")) 
      m_tree_jet->SetBranchAddress("n_jets",            &m_size);
    if (m_tree_jet->FindBranch("jet_4vector")) 
      m_tree_jet->SetBranchAddress("jet_4vector",       &m_jet_lorentzvector);
    if (m_tree_jet->FindBranch("jet_vx")) 
      m_tree_jet->SetBranchAddress("jet_vx",            &m_jet_vx);
    if (m_tree_jet->FindBranch("jet_vy")) 
      m_tree_jet->SetBranchAddress("jet_vy",            &m_jet_vy);
    if (m_tree_jet->FindBranch("jet_vz")) 
      m_tree_jet->SetBranchAddress("jet_vz",            &m_jet_vz);
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
    /*  if (m_tree_jet->FindBranch("jet_btag_BjetProb")) 
        m_tree_jet->SetBranchAddress("jet_btag_BjetProb", &m_jet_btag_BjetProb);
        if (m_tree_jet->FindBranch("jet_btag_SSVHE")) 
        m_tree_jet->SetBranchAddress("jet_btag_SSVHE",    &m_jet_btag_SSVHE);
        if (m_tree_jet->FindBranch("jet_btag_SSVHP")) 
        m_tree_jet->SetBranchAddress("jet_btag_SSVHP",    &m_jet_btag_SSVHP);
        if (m_tree_jet->FindBranch("jet_btag_TCHE")) 
        m_tree_jet->SetBranchAddress("jet_btag_TCHE",     &m_jet_btag_TCHE);  */

    if (m_tree_jet->FindBranch("jet_btag_jetProb")) 
      m_tree_jet->SetBranchAddress("jet_btag_jetProb",  &m_jet_btag_jetProb);
    if (m_tree_jet->FindBranch("jet_btag_TCHP")) 
      m_tree_jet->SetBranchAddress("jet_btag_TCHP",     &m_jet_btag_TCHP);
    if (m_tree_jet->FindBranch("jet_btag_CSV")) 
      m_tree_jet->SetBranchAddress("jet_btag_CSV",      &m_jet_btag_CSV);

    if (m_tree_jet->FindBranch("jet_mcParticleIndex")) 
      m_tree_jet->SetBranchAddress("jet_mcParticleIndex",&m_jet_MCIndex);

    if (m_tree_jet->FindBranch("jet_scaleFactor"))
      m_tree_jet->SetBranchAddress("jet_scaleFactor", &m_scaleFactors.getBackingArray());
  }

  m_tree_met = dynamic_cast<TTree*>(a_file->Get(met_name.c_str()));

  if (!m_tree_met)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_met_lorentzvector = new TClonesArray("TLorentzVector");

  if (m_tree_met->FindBranch("met_4vector"))
    m_tree_met->SetBranchAddress("met_4vector",&m_met_lorentzvector);
}


JetMETExtractor::~JetMETExtractor()
{}


bool JetMETExtractor::isPFJetLoose(const pat::Jet& jet)
{
  // Taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetID

  if (! jet.isPFJet())
    return false;

  // Jet ID works only with uncorrected jet. *EnergyFraction functions take care of that all alone
  bool isValid = jet.neutralHadronEnergyFraction() < 0.99;
  isValid &= jet.neutralEmEnergyFraction() < 0.99;
  isValid &= jet.getPFConstituents().size() > 1;
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
  event.getByLabel(m_tag, jetHandle);
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
  event.getByLabel(m_metTag, metHandle);
  pat::MET MET = (*metHandle).at(0);

  // If we are redoing jet correction, or if the user forces it, recompute TypeI MEt corrections
  if (mCorrectJets || mRedoTypeI) {
    // Raw MET
    edm::Handle<pat::METCollection> rawMetHandle;
    event.getByLabel("patPFMetPFlow", rawMetHandle);

    if (rawMetHandle.isValid()) {
      const pat::MET& rawMet = rawMetHandle->back();
      correctMETWithTypeI(rawMet, MET, p_jets);
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

    if (! isPFJetLoose(rawJet))
      continue;

    JetMETExtractor::writeInfo(event, iSetup, p_jets.at(i), m_size); 

    if (m_MC)
      doMCMatch(p_jets.at(i), event, m_MC, m_size);

    m_size++;
  }

  

  writeInfo(event, iSetup, MET, 0);

  fillTree();
}

void JetMETExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::Jet& part, int index) 
{
  if (index>=m_jets_MAX) return;

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Writing jet #" << index << std::endl;
  std::cout << "Pt: " << part.pt() << "; Px / Pz / Pz / E : " << part.px() << " / " << part.py() << " / " << part.pz() << " / " << part.energy() << std::endl;
  std::cout << "Eta: " << part.eta() << "; Phi : " << part.phi() << std::endl;
#endif

  new((*m_jet_lorentzvector)[index]) TLorentzVector(part.px(),part.py(),part.pz(),part.energy());

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

    //m_jet_btag_jetProb[index]  = part.bDiscriminator("jetProbabilityBJetTags");
    //m_jet_btag_BjetProb[index] = part.bDiscriminator("jetBProbabilityBJetTags");
    //m_jet_btag_SSVHE[index]    = part.bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    //m_jet_btag_SSVHP[index]    = part.bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
    //m_jet_btag_TCHE[index]     = part.bDiscriminator("trackCountingHighEffBJetTags");
    //m_jet_btag_TCHP[index]     = part.bDiscriminator("trackCountingHighPurBJetTags");

    // 2012: See https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
    m_jet_btag_jetProb[index]  = part.bDiscriminator("jetProbabilityBJetTags");
    m_jet_btag_TCHP[index]     = part.bDiscriminator("trackCountingHighPurBJetTags");
    m_jet_btag_CSV[index]      = part.bDiscriminator("combinedSecondaryVertexBJetTags");
  }

  if (m_isMC)
    m_scaleFactors.push_back(m_scaleFactorService->getBTaggingScaleFactor(part.et(), part.eta()));
}

void JetMETExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const pat::MET& part, int index) 
{
  if (index > 1)
    return;

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Writing MET #" << index << std::endl;
  std::cout << "Pt: " << part.pt() << "; Px / Pz / Pz : " << part.px() << " / " << part.py() << " / " << part.pz() << std::endl;
  std::cout << "Eta: " << part.eta() << "; Phi : " << part.phi() << std::endl;
#endif

  new((*m_met_lorentzvector)[index]) TLorentzVector(part.px(),part.py(),part.pz(),part.energy());
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

  m_scaleFactors.clear();

  for (int i=0;i<m_jets_MAX;++i) 
  {
    m_jet_vx[i] = 0.;
    m_jet_vy[i] = 0.;
    m_jet_vz[i] = 0.;

    m_jet_chmult[i] = 0;
    m_jet_chmuEfrac[i] = 0.;
    m_jet_chemEfrac[i] = 0.;
    m_jet_chhadEfrac[i] = 0.;
    m_jet_nemEfrac[i] = 0.;
    m_jet_nhadEfrac[i] = 0.;
    //m_jet_btag_BjetProb[i] = 0.;
    //m_jet_btag_SSVHE[i] = 0.;
    //m_jet_btag_SSVHP[i] = 0.;
    //m_jet_btag_TCHE[i] = 0.;
    m_jet_btag_jetProb[i] = 0.;
    m_jet_btag_TCHP[i] = 0.;
    m_jet_btag_CSV[i] = 0.;

    m_jet_MCIndex[i]    = -1;
  }
  if (m_jet_lorentzvector)
    m_jet_lorentzvector->Clear();

  if (m_met_lorentzvector)
    m_met_lorentzvector->Clear();
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
  const JetCorrector* corrector = JetCorrector::getJetCorrector(mJetCorrectorLabel, iSetup);

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

    double corrections = corrector->correction(jet, iEvent, iSetup);
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

  double factor = 0.;
  double error  = 0.;

  if (fabs(jet.eta()) > 0. && fabs(jet.eta()) <= 0.5) {
    factor = 1.052;
    error = (mJERSign == 1) ? 0.062 : 0.061;
  } else if (fabs(jet.eta()) > 0.5 && fabs(jet.eta()) <= 1.1) {
    factor = 1.057;
    error = (mJERSign == 1) ? 0.056 : 0.055;
  } else if (fabs(jet.eta()) > 1.1 && fabs(jet.eta()) <= 1.7) {
    factor = 1.096;
    error = (mJERSign == 1) ? 0.063 : 0.062;
  } else if (fabs(jet.eta()) > 1.7 && fabs(jet.eta()) <= 2.3) {
    factor = 1.134;
    error = (mJERSign == 1) ? 0.087 : 0.085;
  } else if (fabs(jet.eta()) > 2.3 && fabs(jet.eta()) <= 5.0) {
    factor = 1.288;
    error = (mJERSign == 1) ? 0.155 : 0.153;
  }

  return factor + mJERSign * error;
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

      jet.scaleEnergy(scalefac);
      //rawJet->scaleEnergy(ptscale); 

#if DEBUG
      std::cout << "Corrected pt: " << jet.pt() << std::endl;
#endif
    
      correctedMetPx -= (rawJet->px()*scalefac);
      correctedMetPy -= (rawJet->py()*scalefac);
    }
  }

  
  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);
    
  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));
    

  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorter);
}
    

void JetMETExtractor::correctMETWithTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets) {
  double deltaPx = 0., deltaPy = 0.;
#if DEBUG
    std::cout << "---" << std::endl;
    std::cout << "Computing TypeI correction" << std::endl;
    std::cout << "MET Raw et: " << rawMet.et() << std::endl;
    std::cout << "PAT corrected MET et: " << met.et() << std::endl;
#endif

  // See https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=174324 slide 4
  // and http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h?revision=1.8&view=markup
  for (pat::JetCollection::const_iterator it = jets.begin(); it != jets.end(); ++it) {
    const pat::Jet& jet = *it;

    if (jet.pt() > 10) {

      const pat::Jet* rawJet = jet.userData<pat::Jet>("rawJet");
      const pat::Jet* L1Jet  = jet.userData<pat::Jet>("L1Jet");

      double emEnergyFraction = rawJet->chargedEmEnergyFraction() + rawJet->neutralEmEnergyFraction();
      if (emEnergyFraction > 0.90)
        continue;

      //reco::Candidate::LorentzVector rawJetP4 = rawJet->p4();
      reco::Candidate::LorentzVector L1JetP4  = L1Jet->p4();

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
    }
  }

  double correctedMetPx = rawMet.px() - deltaPx;
  double correctedMetPy = rawMet.py() - deltaPy;
  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);

  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));

#if DEBUG
    std::cout << "Handmade corrected MET et: " << met.et() << std::endl;
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
  event.getByLabel("goodOfflinePrimaryVertices", vertexHandle);
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

#if DEBUG
    std::cout << "---" << std::endl;
    std::cout << "Pt before JES uncertainty: " << jet.pt() << std::endl;
#endif

    jet.scaleEnergy(scaleFactor);

#if DEBUG
    std::cout << "Pt after JES uncertainty: " << jet.pt() << std::endl;
#endif

    correctedMetPx -= (jet.px());
    correctedMetPy -= (jet.py());
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
