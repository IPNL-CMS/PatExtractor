#include "../interface/KVFExtractor.h"
#include "../interface/JECReader.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>

#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#define DEBUG false

using namespace std;

//----------------------------------------------------------------------
// functions from : RecoVertex/KinematicFit/plugins/KineExample.cc

void printOut(const RefCountedKinematicVertex& myVertex)
{
  if (myVertex->vertexIsValid()) {
    cout << "Decay vertex: " << myVertex->position() <<myVertex->chiSquared()<< " "<<myVertex->degreesOfFreedom()<<endl;
  } else cout << "Decay vertex Not valid\n";
}

void printOut(const RefCountedKinematicParticle& myParticle)
{
  cout << "Particle: \n";
  //accessing the reconstructed Bs meson parameters:
  //SK: uncomment if needed  AlgebraicVector7 bs_par = myParticle->currentState().kinematicParameters().vector();

  //and their joint covariance matrix:
  //SK:uncomment if needed  AlgebraicSymMatrix77 bs_er = myParticle->currentState().kinematicParametersError().matrix();
  cout << "Momentum at vertex: " << myParticle->currentState().globalMomentum ()<<endl;
  cout << "Parameters at vertex: " << myParticle->currentState().kinematicParameters().vector()<<endl;
}

void printOut(const RefCountedKinematicTree& myTree)
{
  if (!myTree->isValid()) {
    cout <<"Tree is invalid. Fit failed.\n";
    return;
  }

  //accessing the tree components, move pointer to top
  myTree->movePointerToTheTop();

  //We are now at the top of the decay tree getting the jpsi reconstructed KinematicPartlcle
  RefCountedKinematicParticle jpsi = myTree->currentParticle();
  printOut(jpsi);

  // The jpsi decay vertex
  RefCountedKinematicVertex b_dec_vertex = myTree->currentDecayVertex();
  printOut(b_dec_vertex);

  // Get all the children of Jpsi:
  //In this way, the pointer is not moved
  vector< RefCountedKinematicParticle > jpsi_children = myTree->finalStateParticles();

  for (unsigned int i=0;i< jpsi_children.size();++i) {
    printOut(jpsi_children[i]);
  }

  //Now navigating down the tree , pointer is moved:
  bool child = myTree->movePointerToTheFirstChild();

  if(child) while (myTree->movePointerToTheNextChild()) {
    RefCountedKinematicParticle aChild = myTree->currentParticle();
    printOut(aChild);
  }
}

//-------------------------------------------------------------------------

  KVFExtractor::KVFExtractor(const std::string& name, std::shared_ptr<ScaleFactorService> sf, const edm::ParameterSet& config)
: BaseExtractor(name, sf)
{
  if (! config.exists(name))
    throw edm::Exception(edm::errors::ConfigFileReadError) << "No edm::ParameterSet named " << name << " found";

  const edm::ParameterSet& jetConfig = config.getParameter<edm::ParameterSet>(name);

  m_tag = jetConfig.getParameter<edm::InputTag>("input");

  mCorrectJets = jetConfig.getUntrackedParameter<bool>("redoJetCorrection", false);
  mUseGlobalTagForJEC = jetConfig.getUntrackedParameter<bool>("useGlobalTagForJEC", true);

  if (!mUseGlobalTagForJEC) {
    mJecPayload =  jetConfig.getUntrackedParameter<std::string>("jecPayload");
    mJecJetAlgo =  jetConfig.getUntrackedParameter<std::string>("jecJetAlgo");
    if (mJecPayload.length() > 0) {
      mJecPayload = edm::FileInPath(mJecPayload).fullPath();
    } else {
      std::cout << "WARNING! No JecPayload file found. Use the global tag instead for JEC" << std::endl;
      mUseGlobalTagForJEC = true;
    }
  }

  mTxtCorrector = nullptr;

  if (mCorrectJets)
    mJetCorrectorLabel = jetConfig.getParameter<std::string>("jetCorrectorLabel");
  mDoJER = jetConfig.getUntrackedParameter<bool>("doJER", true);
  mJERSign = 0;
  if (mDoJER) {
    mJERSign = jetConfig.getUntrackedParameter<int>("jerSign", 0);
  }
  mDoLooseJetID = jetConfig.getUntrackedParameter<bool>("doLooseJetID", true);

  mJESSign = jetConfig.getUntrackedParameter<int>("jesSign", 0);

  if (mJERSign != 0 && mJERSign != -1 && mJERSign != 1)
    throw edm::Exception(edm::errors::LogicError) << "jerSign must be 0 for nominal correction, -1 for 1-sigma down correction, or 1 for 1-sigma up correction";

  if (mJESSign != 0 && mJESSign != -1 && mJESSign != 1)
    throw edm::Exception(edm::errors::LogicError) << "jesSign must be 0 for nominal correction, -1 for 1-sigma down correction, or 1 for 1-sigma up correction";

  mJecFilename = jetConfig.getUntrackedParameter<std::string>("jes_uncertainties_file", "");
  if (mJESSign != 0 && mJecFilename.length() > 0) {
    mJecFilename = edm::FileInPath(mJecFilename).fullPath();
  }

#if DEBUG
  std::cout << "##########" << std::endl;
  std::cout << "KVF extractor summary" << std::endl;
  std::cout << "Redo JEC: " << mCorrectJets << ((mCorrectJets) ? ("; " + mJetCorrectorLabel) : "") << std::endl;
  std::cout << "Do JER: " << mDoJER << "; JER sign: " << mJERSign << std::endl;
#endif

  // Set everything to 0
  m_jpsi_jet_lorentzvector     = new TClonesArray("TLorentzVector");
  m_jpsipf_lorentzvector     = new TClonesArray("TLorentzVector");
  m_jpsikvf_lorentzvector     = new TClonesArray("TLorentzVector");
  m_jpsikvf_mu1_lorentzvector = new TClonesArray("TLorentzVector");
  m_jpsikvf_mu2_lorentzvector = new TClonesArray("TLorentzVector");

  reset();

  // Tree definition

  m_tree_jpsi = NULL;
  m_tree_jpsi     = new TTree("jpsi_KVF", "Jpsi info");  

  m_tree_jpsi->Branch("n_jpsi",          &m_jpsi_size,   "n_jpsi/I");

  m_tree_jpsi->Branch("jpsi_indjet",      &m_jpsi_indjet,"jpsi_indjet[n_jpsi]/I");
  m_tree_jpsi->Branch("jpsi_indpf1",      &m_jpsi_indpf1,"jpsi_indpf1[n_jpsi]/I");
  m_tree_jpsi->Branch("jpsi_indpf2",      &m_jpsi_indpf2,"jpsi_indpf2[n_jpsi]/I");

  m_tree_jpsi->Branch("jpsi_jet_4vector",     "TClonesArray",&m_jpsi_jet_lorentzvector, 1000, 0);
  m_tree_jpsi->Branch("jpsipf_4vector",     "TClonesArray",&m_jpsipf_lorentzvector, 1000, 0);

  m_tree_jpsi->Branch("jpsi_4vector",    "TClonesArray",&m_jpsikvf_lorentzvector, 1000, 0);
  m_tree_jpsi->Branch("jpsi_mu1_4vector","TClonesArray",&m_jpsikvf_mu1_lorentzvector, 1000, 0);
  m_tree_jpsi->Branch("jpsi_mu2_4vector","TClonesArray",&m_jpsikvf_mu2_lorentzvector, 1000, 0);
  m_tree_jpsi->Branch("jpsi_vx",	   &m_jpsikvf_vx,	"jpsikvf_vx[n_jpsi]/F");  
  m_tree_jpsi->Branch("jpsi_vy",	   &m_jpsikvf_vy,	"jpsikvf_vy[n_jpsi]/F");  
  m_tree_jpsi->Branch("jpsi_vz",	   &m_jpsikvf_vz,	"jpsikvf_vz[n_jpsi]/F");
  m_tree_jpsi->Branch("jpsi_vtxvalid", &m_jpsikvf_vtxvalid, "jpsikvf_vtxvalid[n_jpsi]/O");
  m_tree_jpsi->Branch("jpsi_vtxchi2",  &m_jpsikvf_vtxchi2,  "jpsikvf_vtxchi2[n_jpsi]/F");
  m_tree_jpsi->Branch("jpsi_ndf",	   &m_jpsikvf_ndf,	"jpsikvf_ndf[n_jpsi]/F");  
  m_tree_jpsi->Branch("jpsi_L3D",	   &m_jpsikvf_L3D,	"jpsikvf_L3D[n_jpsi]/F");  
  m_tree_jpsi->Branch("jpsi_sigmaL3D", &m_jpsikvf_sigmaL3D, "jpsikvf_sigmaL3D[n_jpsi]/F");  
  m_tree_jpsi->Branch("jpsi_L3DoverSigmaL3D", &m_jpsikvf_L3DoverSigmaL3D, "jpsikvf_L3DoverSigmaL3D[n_jpsi]/F");  

}

void KVFExtractor::beginJob() {
  if (!mUseGlobalTagForJEC) {
    mTxtCorrector = makeFactorizedJetCorrectorFromXML(mJecPayload, mJecJetAlgo, m_isMC);
    std::cout << "Using text files for JEC" << std::endl;
  } else {
    std::cout << "Using global tag for JEC" << std::endl;
  }
}


  KVFExtractor::KVFExtractor(const std::string& name, std::shared_ptr<ScaleFactorService> sf, TFile *a_file)
: BaseExtractor(name, sf)
{

  std::cout << "KVFExtractor objet is retrieved" << std::endl;
  m_file = a_file;

  // Tree definition
  m_OK = false;

  m_tree_jpsi = dynamic_cast<TTree*>(a_file->Get(name.c_str()));


  if (m_tree_jpsi) {

    if (m_tree_jpsi->FindBranch("n_jpsi")) 
      m_tree_jpsi->Branch("n_jpsi",          &m_jpsi_size);
    if (m_tree_jpsi->FindBranch("jpsi_indjet")) 
      m_tree_jpsi->Branch("jpsi_indjet",      &m_jpsi_indjet);
    if (m_tree_jpsi->FindBranch("jpsi_indpf1")) 
      m_tree_jpsi->Branch("jpsi_indpf1",      &m_jpsi_indpf1);
    if (m_tree_jpsi->FindBranch("jpsi_indpf2")) 
      m_tree_jpsi->Branch("jpsi_indpf2",      &m_jpsi_indpf2);

    m_jpsi_jet_lorentzvector = new TClonesArray("TLorentzVector");
    m_jpsipf_lorentzvector = new TClonesArray("TLorentzVector");
    m_jpsikvf_mu1_lorentzvector = new TClonesArray("TLorentzVector");
    m_jpsikvf_mu2_lorentzvector = new TClonesArray("TLorentzVector");

    if (m_tree_jpsi->FindBranch("jpsi_jet_4vector")) 
      m_tree_jpsi->Branch("jpsi_jet_4vector",     &m_jpsi_jet_lorentzvector);
    if (m_tree_jpsi->FindBranch("jpsipf_4vector")) 
      m_tree_jpsi->Branch("jpsipf_4vector",     &m_jpsipf_lorentzvector);
    if (m_tree_jpsi->FindBranch("jpsi_4vector")) 
      m_tree_jpsi->Branch("jpsi_4vector",    &m_jpsikvf_lorentzvector);
    if (m_tree_jpsi->FindBranch("jpsi_mu1_4vector")) 
      m_tree_jpsi->Branch("jpsi_mu1_4vector",&m_jpsikvf_mu1_lorentzvector);
    if (m_tree_jpsi->FindBranch("jpsi_mu2_4vector")) 
      m_tree_jpsi->Branch("jpsi_mu2_4vector",&m_jpsikvf_mu2_lorentzvector);
    if (m_tree_jpsi->FindBranch("jpsi_vx")) 
      m_tree_jpsi->Branch("jpsi_vx",	   &m_jpsikvf_vx);  
    if (m_tree_jpsi->FindBranch("jpsi_vy")) 
      m_tree_jpsi->Branch("jpsi_vy",	   &m_jpsikvf_vy);  
    if (m_tree_jpsi->FindBranch("jpsi_vz")) 
      m_tree_jpsi->Branch("jpsi_vz",	   &m_jpsikvf_vz);
    if (m_tree_jpsi->FindBranch("jpsi_vtxvalid")) 
      m_tree_jpsi->Branch("jpsi_vtxvalid", &m_jpsikvf_vtxvalid);
    if (m_tree_jpsi->FindBranch("jpsi_vtxchi2")) 
      m_tree_jpsi->Branch("jpsi_vtxchi2",  &m_jpsikvf_vtxchi2);
    if (m_tree_jpsi->FindBranch("jpsi_ndf")) 
      m_tree_jpsi->Branch("jpsi_ndf",	   &m_jpsikvf_ndf);  
    if (m_tree_jpsi->FindBranch("jpsi_L3D")) 
      m_tree_jpsi->Branch("jpsi_L3D",	   &m_jpsikvf_L3D);  
    if (m_tree_jpsi->FindBranch("jpsi_sigmaL3D")) 
      m_tree_jpsi->Branch("jpsi_sigmaL3D", &m_jpsikvf_sigmaL3D);  
    if (m_tree_jpsi->FindBranch("jpsi_L3DoverSigmaL3D")) 
      m_tree_jpsi->Branch("jpsi_L3DoverSigmaL3D", &m_jpsikvf_L3DoverSigmaL3D);  
  }

  m_OK = true;

}


KVFExtractor::~KVFExtractor()
{
  delete mTxtCorrector;
}


bool KVFExtractor::isPFJetLoose(const pat::Jet& jet)
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

void KVFExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* m_MC) 
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

  //do Jets MET resolution corrections
  if (m_MC && mDoJER)
    correctJetsResolution(p_jets);

  // JES systematics
  if (mJESSign != 0)
    doJESSystematics(p_jets);

  int nJpsi = 0; 

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

    // Reconstruct the J/psi
    //----------------------

    const std::vector<reco::PFCandidatePtr>& PFpart = p_jets.at(i).getPFConstituents();
    unsigned int npfs = PFpart.size(); 

    ParticleMass muon_mass  = 0.1056583;
    float        muon_sigma = 0.0000001;

    // ParticleMass jpsi_mass  = 3.09687;
    // float        jpsi_sigma = 0.00009;

    // To transform Track to TransientTrack, first need to get the builder:
    edm::ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);  

    for(unsigned int j = 0; j < npfs; ++j) {

      if (abs(PFpart[j]->pdgId()) != 13) continue;
      if (PFpart[j]->pt() < 4) continue;
      if (!PFpart[j]->trackRef()) continue;

      for(unsigned int k = j+1; k<npfs; ++k) {

        // Both PF particle should be a muon and with OS
        if (abs(PFpart[k]->pdgId()) != 13) continue;
        if (PFpart[k]->pt() < 4) continue;
        if (!PFpart[k]->trackRef()) continue;
        if (PFpart[j]->charge() + PFpart[k]->charge() != 0) continue;

        double eJpsi  = PFpart[j]->energy()+PFpart[k]->energy();
        double pxJpsi = PFpart[j]->px()+PFpart[k]->px();
        double pyJpsi = PFpart[j]->py()+PFpart[k]->py();
        double pzJpsi = PFpart[j]->pz()+PFpart[k]->pz();
        double mJpsi = pow(eJpsi,2)-pow(pxJpsi,2)-pow(pyJpsi,2)-pow(pzJpsi,2);
        if (mJpsi > 0.) mJpsi = sqrt(mJpsi);
        else mJpsi = 0.;

        if (mJpsi >= 2.8 && mJpsi < 3.4) {
          if (nJpsi > 0) std::cout << "Warning: another Jpsi found with this PF particle" << std::endl;
          ++nJpsi;

          // Fill tree info for links between Jpsi and jet
          m_jpsi_indjet[nJpsi-1] = i;
          new((*m_jpsi_jet_lorentzvector)[nJpsi-1]) TLorentzVector((p_jets.at(i)).px(),(p_jets.at(i)).py(),(p_jets.at(i)).pz(),(p_jets.at(i)).energy());
          // Fill tree info for links between Jpsi and PF particles
          m_jpsi_indpf1[nJpsi-1] = j;
          m_jpsi_indpf2[nJpsi-1] = k;
          // save the Jpsi from 2 PF muons (which can be different from Jpsi from 2 tracks)
          new((*m_jpsipf_lorentzvector)[nJpsi-1]) TLorentzVector(pxJpsi, pyJpsi, pzJpsi, eJpsi);

          // Make the Transient tracks
          reco::TransientTrack tr1 = (*theB).build(PFpart[j]->trackRef());
          reco::TransientTrack tr2 = (*theB).build(PFpart[k]->trackRef());

          // A complicated Kalman fit :
          //---------------------------

          //Creating a KinematicParticleFactory
          KinematicParticleFactoryFromTransientTrack pFactory;

          //initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered
          float chi = 0.;
          float ndf = 0.;

          //making particles
          std::vector<RefCountedKinematicParticle> muonParticles;
          muonParticles.push_back(pFactory.particle (tr1,muon_mass,chi,ndf,muon_sigma));
          muonParticles.push_back(pFactory.particle (tr2,muon_mass,chi,ndf,muon_sigma));

          /* Example of a simple vertex fit, without other constraints
           * The reconstructed decay tree is a result of the kinematic fit
           * The KinematicParticleVertexFitter fits the final state particles to their vertex and
           * reconstructs the decayed state
           */

          // creating the vertex fitter
          KinematicParticleVertexFitter fitter;
          // reconstructing a J/Psi decay
          RefCountedKinematicTree vertexFitTree = fitter.fit(muonParticles);

          if (!vertexFitTree->isValid()) {
            cout <<" vertexTree is invalid. Fit failed.\n";

            // Need to fill empty quantities : OK for tables, but need to create empty TLorentzVector
            new((*m_jpsikvf_lorentzvector)[nJpsi-1]) TLorentzVector(0,0,0,0);
            new((*m_jpsikvf_mu1_lorentzvector)[nJpsi-1]) TLorentzVector(0,0,0,0);
            new((*m_jpsikvf_mu2_lorentzvector)[nJpsi-1]) TLorentzVector(0,0,0,0);

            continue;

          } else {
            //accessing the tree components, move pointer to top
            vertexFitTree->movePointerToTheTop();

            //We are now at the top of the decay tree getting the jpsi reconstructed KinematicPartlcle
            RefCountedKinematicParticle jpsi1 = vertexFitTree->currentParticle();
            AlgebraicVector7 par0 = jpsi1->currentState().kinematicParameters().vector();
            double e0 = jpsi1->currentState().kinematicParameters().energy();
            new((*m_jpsikvf_lorentzvector)[nJpsi-1]) TLorentzVector(par0(3),par0(4),par0(5),e0);

            RefCountedKinematicVertex jpsi1_vertex = vertexFitTree->currentDecayVertex();
            if ( jpsi1_vertex->vertexIsValid()) {

              m_jpsikvf_vx[nJpsi-1] = jpsi1_vertex->position().x();
              m_jpsikvf_vy[nJpsi-1] = jpsi1_vertex->position().y();
              m_jpsikvf_vz[nJpsi-1] = jpsi1_vertex->position().z();
              m_jpsikvf_vtxvalid[nJpsi-1] = true;
              m_jpsikvf_vtxchi2[nJpsi-1] = jpsi1_vertex->chiSquared();
              m_jpsikvf_ndf[nJpsi-1] = jpsi1_vertex->degreesOfFreedom();

              // Compute the distance between the PV and the Jpsi vertex :
              //----------------------------------------------------------

              edm::Handle<reco::VertexCollection>  vtxHandle;
              edm::InputTag tagVtx("offlinePrimaryVertices");
              event.getByLabel(tagVtx, vtxHandle);
              const reco::VertexCollection vtx = *(vtxHandle.product());

              GlobalPoint svPos    = jpsi1_vertex->position();
              GlobalError svPosErr = jpsi1_vertex->error();

              // If the  PV does not exist, compute distance wrt to the detector center (0,0,0)

              double sigmax = 0.;
              double sigmay = 0.;
              double sigmaz = 0.;
              if (vtx.size() > 0) {
                sigmax = sqrt( vtx[0].xError()*vtx[0].xError() + svPosErr.cxx()*svPosErr.cxx() );
                sigmay = sqrt( vtx[0].yError()*vtx[0].yError() + svPosErr.cyy()*svPosErr.cyy() );
                sigmaz = sqrt( vtx[0].zError()*vtx[0].zError() + svPosErr.czz()*svPosErr.czz() );
              }  else {
                sigmax = sqrt( svPosErr.cxx()*svPosErr.cxx() );
                sigmay = sqrt( svPosErr.cyy()*svPosErr.cyy() );
                sigmaz = sqrt( svPosErr.czz()*svPosErr.czz() );
              }

              double px  = ((TLorentzVector*)m_jpsikvf_lorentzvector->At(nJpsi-1))->Px();
              double py  = ((TLorentzVector*)m_jpsikvf_lorentzvector->At(nJpsi-1))->Py();
              double pz  = ((TLorentzVector*)m_jpsikvf_lorentzvector->At(nJpsi-1))->Pz();
              double nrj = ((TLorentzVector*)m_jpsikvf_lorentzvector->At(nJpsi-1))->E();
              double m = sqrt( nrj*nrj - px*px - py*py - pz*pz );

              double interx = pow( (px/m)/sigmax, 2.);
              double intery = pow( (py/m)/sigmay, 2.);
              double interz = pow( (pz/m)/sigmaz, 2.);

              m_jpsikvf_sigmaL3D[nJpsi-1] = pow( interx + intery + interz , -0.5);
              //std::cout << "sigmaL3D = " << m_jpsikvf_sigmaL3D[nJpsi-1] << std::endl;

              double part1 = 0.;
              double part2 = 0.;
              double part3 = 0.;
              if (vtx.size() >0) {
                part1 = (px/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmax,2.)*( svPos.x() - vtx[0].x());
                part2 = (py/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmay,2.)*( svPos.y() - vtx[0].y());
                part3 = (pz/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmaz,2.)*( svPos.z() - vtx[0].z());
              }  else {
                part1 = (px/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmax,2.)*svPos.x();
                part2 = (py/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmay,2.)*svPos.y();
                part3 = (pz/m)*pow(m_jpsikvf_sigmaL3D[nJpsi-1]/sigmaz,2.)*svPos.z();
              }

              m_jpsikvf_L3D[nJpsi-1] = fabs(part1 + part2 + part3);
              //std::cout << "L3D = " << m_jpsikvf_L3D[nJpsi-1] << std::endl;  

              m_jpsikvf_L3DoverSigmaL3D[nJpsi-1] = m_jpsikvf_L3D[nJpsi-1]/m_jpsikvf_sigmaL3D[nJpsi-1];
              //std::cout << "(L/sigma)3D = " << m_jpsikvf_L3DoverSigmaL3D[nJpsi-1] << std::endl;    

              if (!vtxHandle.isValid()) {
                std::cout << "KVFExtractor::writeInfo(): vtxHandle is not valid..." << std::endl;
              }

              vector< RefCountedKinematicParticle > jpsi1_children = vertexFitTree->finalStateParticles();
              if (jpsi1_children.size() != 2) {
                cout << " Warning Jpsi1 children size not equal to 2..." << endl;
              } else {
                // Order is : x,y,z,px,py,pz,m
                AlgebraicVector7 par1 = jpsi1_children[0]->currentState().kinematicParameters().vector();
                double e1 = jpsi1_children[0]->currentState().kinematicParameters().energy();
                new((*m_jpsikvf_mu1_lorentzvector)[nJpsi-1]) TLorentzVector(par1(3),par1(4),par1(5),e1);

                AlgebraicVector7 par2 = jpsi1_children[1]->currentState().kinematicParameters().vector();
                double e2 = jpsi1_children[1]->currentState().kinematicParameters().energy();
                new((*m_jpsikvf_mu2_lorentzvector)[nJpsi-1]) TLorentzVector(par2(3),par2(4),par2(5),e2);
              }

            } else cout << "Decay vertex Not valid\n";

          } // vertex is not valid
        } // mass condition

      } // end 2nd PF loop
    } // end 1st PF loop

    m_jpsi_size = nJpsi;

    // end of J/psi stuff

    m_size++;
  }  // jet loop

  fillTree();
}

//
// Method getting the info from an input file
//

void KVFExtractor::getInfo(int ievt) 
{
  if (m_tree_jpsi)
    m_tree_jpsi->GetEntry(ievt);

}


// Method initializing everything (to do for each event)

void KVFExtractor::reset()
{
  m_size = 0;


  // Jpsi tree
  m_jpsi_jet_lorentzvector->Clear();
  m_jpsipf_lorentzvector->Clear();
  m_jpsikvf_lorentzvector->Clear();
  m_jpsikvf_mu1_lorentzvector->Clear();
  m_jpsikvf_mu2_lorentzvector->Clear();

  m_jpsi_size = 0;

  for (int i=0;i<m_jpsi_MAX;++i) {
    m_jpsi_indjet[i] = 0;
    m_jpsi_indpf1[i] = 0;
    m_jpsi_indpf2[i] = 0;

    m_jpsikvf_vx[i] = 0.;
    m_jpsikvf_vy[i] = 0.;
    m_jpsikvf_vz[i] = 0.;
    m_jpsikvf_vtxvalid[i] = false;
    m_jpsikvf_vtxchi2[i] = 0.;
    m_jpsikvf_ndf[i] = 0.;

    m_jpsikvf_L3D[i] = 0.;
    m_jpsikvf_sigmaL3D[i] = 0.;
    m_jpsikvf_L3DoverSigmaL3D[i] = 0.;
  }
}


void KVFExtractor::fillTree()
{
  if (m_tree_jpsi)
    m_tree_jpsi->Fill(); 

}

void KVFExtractor::correctJets(pat::JetCollection& jets, const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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
  iEvent.getByLabel("goodOfflinePrimaryVertices", vertexHandle);
  reco::VertexCollection vertices = *vertexHandle;

  edm::Handle<double> rhos;
  iEvent.getByLabel(edm::InputTag("kt6PFJets", "rho", "RECO"), rhos);
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

double KVFExtractor::getResCorrFactor(const pat::Jet& jet) {

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

void KVFExtractor::correctJetsResolution(pat::JetCollection& jets) {

#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "Doing jet resolution smearing" << std::endl;
#endif

  // Correct jets
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    pat::Jet& jet = *it;

    if (jet.pt() > 10) {

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

      jet.scaleEnergy(scalefac);

#if DEBUG
      std::cout << "Corrected pt: " << jet.pt() << std::endl;
#endif

    }
  }

  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorter);
}


//from JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi, check updates
//updated to MET systematic shift corrections to 2012 ABCD ReReco data + new Summer'13 JEC
double KVFExtractor::getSysShifCorrFactorX(const int Nvtx){
  if (m_isMC) return -(+1.62861e-01 - 2.38517e-02*Nvtx);
  else        return -(+4.83642e-02 + 2.48870e-01*Nvtx);
}

double KVFExtractor::getSysShifCorrFactorY(const int Nvtx){
  if (m_isMC) return -(+3.60860e-01 - 1.30335e-01*Nvtx);
  else        return -(-1.50135e-01 - 8.27917e-02*Nvtx);
}


void KVFExtractor::doJESSystematics(pat::JetCollection& jets) {
#if DEBUG
  std::cout << "---" << std::endl;
  std::cout << "JES systematic" << std::endl;
#endif

  // Correct jets
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    pat::Jet& jet = *it;

    jecUncertainty->setJetEta(jet.eta());
    jecUncertainty->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt

    double uncertainty = (mJESSign == 1) ? fabs(jecUncertainty->getUncertainty(true)) : fabs(jecUncertainty->getUncertainty(false));
    double signedCorrection = mJESSign * uncertainty;

    double scaleFactor = 1. + signedCorrection;

#if DEBUG
    std::cout << "---" << std::endl;
    std::cout << "Pt before JES uncertainty: " << jet.pt() << std::endl;
#endif

    jet.scaleEnergy(scaleFactor);

#if DEBUG
    std::cout << "Pt after JES uncertainty: " << jet.pt() << std::endl;
#endif

  }

  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorter);
}

void KVFExtractor::extractRawJets(pat::JetCollection& jets) {

  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it) {
    pat::Jet& jet = *it;

    const pat::Jet rawJet = jet.correctedJet("Uncorrected");
    jet.addUserData("rawJet", rawJet, true);
    const pat::Jet L1Jet  = jet.correctedJet("L1FastJet");
    jet.addUserData("L1Jet", L1Jet, true); // Embed L1 corrected jet for TypeI correction
  }

}

