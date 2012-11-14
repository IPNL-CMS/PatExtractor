#include "../interface/JetMETExtractor.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#define DEBUG false

JetMETExtractor::JetMETExtractor(const std::string& name, const std::string& met_name, const edm::InputTag& tag, const edm::InputTag& metTag,
    bool doJetTree, bool doMETTree, bool correctJets, const std::string& jetCorrectorLabel, bool redoTypeI)
: BaseExtractor(name)
{
  m_tag = tag;
  m_metTag = metTag;

  mCorrectJets = correctJets;
  mJetCorrectorLabel = jetCorrectorLabel;
  mRedoTypeI = redoTypeI;

  // Set everything to 0
  m_jet_lorentzvector = new TClonesArray("TLorentzVector");
  m_met_lorentzvector = new TClonesArray("TLorentzVector");

  reset();

  // Tree definition

  m_tree_jet = NULL;
  if (doJetTree)
  {
    m_tree_jet     = new TTree(name.c_str(), "PAT PF jet info");  
    m_tree_jet->Branch("n_jets",  &m_size,   "n_jets/I");  
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
  }

  m_tree_met = NULL;
  if (doMETTree) {
    m_tree_met      = new TTree(met_name.c_str(), "PAT PF MET info");  
    m_tree_met->Branch("met_4vector","TClonesArray",&m_met_lorentzvector, 1000, 0);
  }
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

  if (mCorrectJets) {
    correctJets(p_jets, event, iSetup);
  } else {
    extractRawJets(p_jets);
  }

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

    JetMETExtractor::writeInfo(p_jets.at(i), m_size); 

    if (m_MC)
      doMCMatch(p_jets.at(i), event, m_MC, m_size);

    m_size++;
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

  writeInfo(MET, 0);

  fillTree();
}

void JetMETExtractor::writeInfo(const pat::Jet& part, int index) 
{
  if (index>=m_jets_MAX) return;

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
}

void JetMETExtractor::writeInfo(const pat::MET& part, int index) 
{
  if (index > 1)
    return;

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

void JetMETExtractor::correctMETWithTypeI(const pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets) {
  double deltaPx = 0., deltaPy = 0.;
#if DEBUG
    std::cout << "---" << std::endl;
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
