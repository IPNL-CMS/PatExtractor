#include "../interface/PFpartExtractor.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
// #include "RecoVertex/KinematicFitPrimitives/interface/"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

using namespace std;

//----------------------------------------------------------------------
// functions from : RecoVertex/KinematicFit/plugins/KineExample.cc

void printout(const RefCountedKinematicVertex& myVertex)
{
  if (myVertex->vertexIsValid()) {
    cout << "Decay vertex: " << myVertex->position() <<myVertex->chiSquared()<< " "<<myVertex->degreesOfFreedom()<<endl;
  } else cout << "Decay vertex Not valid\n";
}

void printout(const RefCountedKinematicParticle& myParticle)
{
  cout << "Particle: \n";
  //accessing the reconstructed Bs meson parameters:
  //SK: uncomment if needed  AlgebraicVector7 bs_par = myParticle->currentState().kinematicParameters().vector();

  //and their joint covariance matrix:
  //SK:uncomment if needed  AlgebraicSymMatrix77 bs_er = myParticle->currentState().kinematicParametersError().matrix();
  cout << "Momentum at vertex: " << myParticle->currentState().globalMomentum ()<<endl;
  cout << "Parameters at vertex: " << myParticle->currentState().kinematicParameters().vector()<<endl;
}

void printout(const RefCountedKinematicTree& myTree)
{
  if (!myTree->isValid()) {
    cout <<"Tree is invalid. Fit failed.\n";
    return;
  }

  //accessing the tree components, move pointer to top
  myTree->movePointerToTheTop();

  //We are now at the top of the decay tree getting the jpsi reconstructed KinematicPartlcle
  RefCountedKinematicParticle jpsi = myTree->currentParticle();
  printout(jpsi);

  // The jpsi decay vertex
  RefCountedKinematicVertex b_dec_vertex = myTree->currentDecayVertex();
  printout(b_dec_vertex);

  // Get all the children of Jpsi:
  //In this way, the pointer is not moved
  vector< RefCountedKinematicParticle > jpsi_children = myTree->finalStateParticles();

  for (unsigned int i=0;i< jpsi_children.size();++i) {
    printout(jpsi_children[i]);
  }

  //Now navigating down the tree , pointer is moved:
  bool child = myTree->movePointerToTheFirstChild();

  if(child) while (myTree->movePointerToTheNextChild()) {
    RefCountedKinematicParticle aChild = myTree->currentParticle();
    printout(aChild);
  }
}

//-------------------------------------------------------------------------

  PFpartExtractor::PFpartExtractor(const std::string& name, const edm::InputTag& tag, bool doTree)
: BaseExtractor(name)
{
  m_tag = tag;

  // Set everything to 0
  m_OK = false;

  m_pf_lorentzvector = new TClonesArray("TLorentzVector");

  m_jpsimu_lorentzvector     = new TClonesArray("TLorentzVector");

  m_jpsiraw_lorentzvector     = new TClonesArray("TLorentzVector");
  m_jpsiraw_mu1_lorentzvector = new TClonesArray("TLorentzVector");
  m_jpsiraw_mu2_lorentzvector = new TClonesArray("TLorentzVector");

  reset();

  // Tree definition

  if (doTree)
  {
    m_OK = true;
    m_tree_pfpart         = new TTree(name.c_str(), "PF particles info");     
    m_tree_pfpart->Branch("n_pf",      &m_pf_size,   "n_pf/I");  
    m_tree_pfpart->Branch("pf_4vector","TClonesArray",&m_pf_lorentzvector, 1000, 0);
    m_tree_pfpart->Branch("pf_vx",     &m_pf_vx,     "pf_vx[n_pf]/F");  
    m_tree_pfpart->Branch("pf_vy",     &m_pf_vy,     "pf_vy[n_pf]/F");  
    m_tree_pfpart->Branch("pf_vz",     &m_pf_vz,     "pf_vz[n_pf]/F");  
    m_tree_pfpart->Branch("pf_charge", &m_pf_charge, "pf_charge[n_pf]/I");
    m_tree_pfpart->Branch("pf_pdgid",  &m_pf_pdgid,  "pf_charge[n_pf]/I");
    m_tree_pfpart->Branch("pf_trkLayer", &m_pf_trkLayer,  "pf_trkLayer[n_pf]/I");
    m_tree_pfpart->Branch("pf_pixLayer", &m_pf_pixLayer,  "pf_pixLayer[n_pf]/I");
    m_tree_pfpart->Branch("pf_trknormChi2", &m_pf_trknormChi2,  "pf_trknormChi2[n_pf]/F");
    m_tree_pfpart->Branch("pf_numberOfChambers", &m_pf_numberOfChambers,  "pf_numberOfChambers[n_pf]/I");
    m_tree_pfpart->Branch("pf_numberOfMatchedStations", &m_pf_numberOfMatchedStations,  "pf_numberOfMatchedStations[n_pf]/I");
    m_tree_pfpart->Branch("pf_isGlobalMuon", &m_pf_isGlobalMuon,  "pf_isGlobalMuon[n_pf]/O");
    m_tree_pfpart->Branch("pf_isTrackerMuon", &m_pf_isTrackerMuon,  "pf_isTrackerMuon[n_pf]/O");
    m_tree_pfpart->Branch("pf_isStandAloneMuon", &m_pf_isStandAloneMuon,  "pf_isStandAloneMuon[n_pf]/O");
    m_tree_pfpart->Branch("pf_isCaloMuon", &m_pf_isCaloMuon,  "pf_isCaloMuon[n_pf]/O");
    m_tree_pfpart->Branch("pf_isPFMuon", &m_pf_isPFMuon,  "pf_isPFMuon[n_pf]/O");
    m_tree_pfpart->Branch("pf_isRPCMuon", &m_pf_isRPCMuon,  "pf_isRPCMuon[n_pf]/O");

    m_tree_pfpart->Branch("n_jpsi",          &m_jpsi_size,   "n_jpsi/I");

    m_tree_pfpart->Branch("jpsi_indpf1",      &m_jpsi_indpf1,"jpsi_indpf1[n_jpsi]/I");
    m_tree_pfpart->Branch("jpsi_indpf2",      &m_jpsi_indpf2,"jpsi_indpf2[n_jpsi]/I");

    m_tree_pfpart->Branch("jpsimu_4vector",     "TClonesArray",&m_jpsimu_lorentzvector, 1000, 0);

    m_tree_pfpart->Branch("jpsi_4vector",    "TClonesArray",&m_jpsiraw_lorentzvector, 1000, 0);
    m_tree_pfpart->Branch("jpsi_mu1_4vector","TClonesArray",&m_jpsiraw_mu1_lorentzvector, 1000, 0);
    m_tree_pfpart->Branch("jpsi_mu2_4vector","TClonesArray",&m_jpsiraw_mu2_lorentzvector, 1000, 0);
    m_tree_pfpart->Branch("jpsi_vx",	   &m_jpsiraw_vx,	"jpsiraw_vx[n_jpsi]/F");  
    m_tree_pfpart->Branch("jpsi_vy",	   &m_jpsiraw_vy,	"jpsiraw_vy[n_jpsi]/F");  
    m_tree_pfpart->Branch("jpsi_vz",	   &m_jpsiraw_vz,	"jpsiraw_vz[n_jpsi]/F");
    m_tree_pfpart->Branch("jpsi_vtxvalid", &m_jpsiraw_vtxvalid, "jpsiraw_vtxvalid[n_jpsi]/O");
    m_tree_pfpart->Branch("jpsi_vtxchi2",  &m_jpsiraw_vtxchi2,  "jpsiraw_vtxchi2[n_jpsi]/F");
    m_tree_pfpart->Branch("jpsi_ndf",	   &m_jpsiraw_ndf,	"jpsiraw_ndf[n_jpsi]/F");  
    m_tree_pfpart->Branch("jpsi_L3D",	   &m_jpsiraw_L3D,	"jpsiraw_L3D[n_jpsi]/F");  
    m_tree_pfpart->Branch("jpsi_sigmaL3D", &m_jpsiraw_sigmaL3D, "jpsiraw_sigmaL3D[n_jpsi]/F");  
    m_tree_pfpart->Branch("jpsi_L3DoverSigmaL3D", &m_jpsiraw_L3DoverSigmaL3D, "jpsiraw_L3DoverSigmaL3D[n_jpsi]/F");  

  }
}

  PFpartExtractor::PFpartExtractor(const std::string& name, TFile *a_file)
:BaseExtractor(name)
{
  std::cout << "PFpartExtractor objet is retrieved" << std::endl;

  // Tree definition
  m_OK = false;

  m_tree_pfpart = dynamic_cast<TTree*>(a_file->Get(m_name.c_str()));

  if (!m_tree_pfpart)
  {
    std::cout << "This tree doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;

  m_pf_lorentzvector = new TClonesArray("TLorentzVector");

  m_tree_pfpart->SetBranchAddress("n_pf",     &m_pf_size);
  m_tree_pfpart->SetBranchAddress("pf_4vector",&m_pf_lorentzvector);
  m_tree_pfpart->SetBranchAddress("pf_vx",     &m_pf_vx);
  m_tree_pfpart->SetBranchAddress("pf_vy",     &m_pf_vy);
  m_tree_pfpart->SetBranchAddress("pf_vz",     &m_pf_vz);
  m_tree_pfpart->SetBranchAddress("pf_charge", &m_pf_charge);
  m_tree_pfpart->SetBranchAddress("pf_pdgid",  &m_pf_pdgid);

  // Need to do the same for the Jpsi tree

}

PFpartExtractor::~PFpartExtractor()
{}



//
// Method filling the main particle tree
//

// Dummy pure virtual function
void PFpartExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, const reco::PFCandidate& part, int index) 
{

}

void PFpartExtractor::writeInfo(const edm::Event& event, const edm::EventSetup& iSetup, MCExtractor* mcExtractor)
{
  edm::Handle<reco::PFCandidateCollection>  pfHandle;
  event.getByLabel(m_tag, pfHandle);
  reco::PFCandidateCollection pfs = *pfHandle;
  if ( !pfHandle.isValid() ) {
    std::cout << "PFpartExtractor::writeInfo(): pfHandle is not valid..." << std::endl;
    return;
  }

  // Select good PF muons :
  //-----------------------

  reset();

  m_size = 0;

  std::vector<const reco::PFCandidate*> myPFparts;
  std::vector<int> indpfs;

  for ( unsigned int i = 0; i < pfs.size(); ++i ) {

    if ( abs(pfs[i].pdgId()) != 13  ) continue;
    if ( pfs[i].pt()         <   4. ) continue;

    myPFparts.push_back(&pfs[i]);

    if (m_pf_size>=m_pfpart_MAX) continue;

    new((*m_pf_lorentzvector)[m_pf_size]) TLorentzVector(pfs[i].px(),pfs[i].py(),pfs[i].pz(),pfs[i].energy());
    if ( !pfs[i].trackRef() ) {
      m_pf_vx[m_pf_size]              = 0;
      m_pf_vy[m_pf_size]              = 0;
      m_pf_vz[m_pf_size]              = 0;
      m_pf_trkLayer[m_pf_size]        = 0;
      m_pf_pixLayer[m_pf_size]        = 0;
    } else {
      m_pf_vx[m_pf_size]              = pfs[i].trackRef()->vx();
      m_pf_vy[m_pf_size]              = pfs[i].trackRef()->vy();
      m_pf_vz[m_pf_size]              = pfs[i].trackRef()->vz();
      m_pf_trkLayer[m_pf_size]        = pfs[i].trackRef()->hitPattern().trackerLayersWithMeasurement();
      m_pf_pixLayer[m_pf_size]        = pfs[i].trackRef()->hitPattern().pixelLayersWithMeasurement();
      m_pf_trknormChi2[m_pf_size]     = pfs[i].trackRef()->normalizedChi2();
    }

    if ( !pfs[i].muonRef() ) {
      cout << "PF muonRef is invalid ... " << endl;
    } else {
      //      m_pf_numberOfChambers[m_pf_size]        = pfs[i].muonRef()->numberOfChambers();
      //      m_pf_numberOfMatchedStations[m_pf_size] = pfs[i].muonRef()->numberOfMatchedStations();
      //      m_pf_isGlobalMuon[m_pf_size]     =  pfs[i].muonRef()->isGlobalMuon();
      //      m_pf_isTrackerMuon[m_pf_size]    =  pfs[i].muonRef()->isTrackerMuon();
      //      m_pf_isStandAloneMuon[m_pf_size] =  pfs[i].muonRef()->isStandAloneMuon();
      //      m_pf_isCaloMuon[m_pf_size]       =  pfs[i].muonRef()->isCaloMuon();
      //      m_pf_isPFMuon[m_pf_size]         =  pfs[i].muonRef()->isPFMuon();
      //      m_pf_isRPCMuon[m_pf_size]        =  pfs[i].muonRef()->isRPCMuon();
    }

    m_pf_charge[m_pf_size]          = pfs[i].charge();
    m_pf_pdgid[m_pf_size]           = pfs[i].pdgId();


    indpfs.push_back(m_pf_size);
    ++m_pf_size;
  }

  // Reconstruct the J/psi
  //----------------------

  ParticleMass muon_mass  = 0.1056583;
  float        muon_sigma = 0.0000001;

  //ParticleMass jpsi_mass  = 3.09687;
  //float        jpsi_sigma = 0.00009;

  AddFourMomenta addp4;
  std::vector<reco::CompositeCandidate> allJpsiCands;
  std::vector<const reco::PFCandidate*> pfpart1;
  std::vector<const reco::PFCandidate*> pfpart2;

  for(unsigned int j=0; j<myPFparts.size(); ++j ) {

    bool foundAjpsi = false;

    for(unsigned int k=j+1; k<myPFparts.size(); ++k ) {

      // Both PF particle should be a muon and with OS
      if ( abs(myPFparts[j]->pdgId()) != 13 && abs(myPFparts[k]->pdgId()) != 13 ) continue;
      if ( myPFparts[j]->charge()+myPFparts[k]->charge() != 0                   ) continue;
      if ( !myPFparts[j]->trackRef() || !myPFparts[k]->trackRef()               ) continue;

      reco::CompositeCandidate jpsicand;
      jpsicand.addDaughter(*myPFparts[j],"mu1");
      jpsicand.addDaughter(*myPFparts[k],"mu2");
      addp4.set(jpsicand);

      if ( jpsicand.mass() >= 2.8 && jpsicand.mass() < 3.4 ) {
        if ( foundAjpsi ) std::cout << "Warning: another Jpsi found with this PF particle" << std::endl;
        foundAjpsi=true;
        allJpsiCands.push_back(jpsicand);
        pfpart1.push_back(myPFparts[j]);
        pfpart2.push_back(myPFparts[k]);

        // Fill tree info for links between Jpsi and PF particles
        m_jpsi_indpf1[allJpsiCands.size()-1] = indpfs[j];
        m_jpsi_indpf2[allJpsiCands.size()-1] = indpfs[k];
      }

    } // end 2nd loop
  } // end 1st loop

  //  std::cout << "Number of J/psi found = " << allJpsiCands.size() << std::endl;

  // To transform Track to TransientTrack, first need to get the builder:

  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);  

  m_jpsi_size = allJpsiCands.size();

  for ( unsigned int i=0; i<allJpsiCands.size(); ++i) {

    // save the Jpsi from 2 PF muons (which can be different from Jpsi from 2 tracks)
    new((*m_jpsimu_lorentzvector)[i]) TLorentzVector(allJpsiCands[i].px(),allJpsiCands[i].py(),allJpsiCands[i].pz(),allJpsiCands[i].energy());

    // Make the Transient tracks
    reco::TransientTrack tr1 = (*theB).build(pfpart1[i]->trackRef());
    reco::TransientTrack tr2 = (*theB).build(pfpart2[i]->trackRef());

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
      new((*m_jpsiraw_lorentzvector)[i]) TLorentzVector(0,0,0,0);
      new((*m_jpsiraw_mu1_lorentzvector)[i]) TLorentzVector(0,0,0,0);
      new((*m_jpsiraw_mu2_lorentzvector)[i]) TLorentzVector(0,0,0,0);

      continue;

    } else {
      //accessing the tree components, move pointer to top
      vertexFitTree->movePointerToTheTop();

      //We are now at the top of the decay tree getting the jpsi reconstructed KinematicPartlcle
      RefCountedKinematicParticle jpsi1 = vertexFitTree->currentParticle();
      AlgebraicVector7 par0 = jpsi1->currentState().kinematicParameters().vector();
      double e0 = jpsi1->currentState().kinematicParameters().energy();
      new((*m_jpsiraw_lorentzvector)[i]) TLorentzVector(par0(3),par0(4),par0(5),e0);

      RefCountedKinematicVertex jpsi1_vertex = vertexFitTree->currentDecayVertex();
      if ( jpsi1_vertex->vertexIsValid()) {

        m_jpsiraw_vx[i] = jpsi1_vertex->position().x();
        m_jpsiraw_vy[i] = jpsi1_vertex->position().y();
        m_jpsiraw_vz[i] = jpsi1_vertex->position().z();
        m_jpsiraw_vtxvalid[i] = true;
        m_jpsiraw_vtxchi2[i]  = jpsi1_vertex->chiSquared();
        m_jpsiraw_ndf[i]      = jpsi1_vertex->degreesOfFreedom();

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
        if (vtx.size() >0) {
          sigmax = sqrt( vtx[0].xError()*vtx[0].xError() + svPosErr.cxx()*svPosErr.cxx() );
          sigmay = sqrt( vtx[0].yError()*vtx[0].yError() + svPosErr.cyy()*svPosErr.cyy() );
          sigmaz = sqrt( vtx[0].zError()*vtx[0].zError() + svPosErr.czz()*svPosErr.czz() );
        }  else {
          sigmax = sqrt( svPosErr.cxx()*svPosErr.cxx() );
          sigmay = sqrt( svPosErr.cyy()*svPosErr.cyy() );
          sigmaz = sqrt( svPosErr.czz()*svPosErr.czz() );
        }

        double px  = ((TLorentzVector*)m_jpsiraw_lorentzvector->At(i))->Px();
        double py  = ((TLorentzVector*)m_jpsiraw_lorentzvector->At(i))->Py();
        double pz  = ((TLorentzVector*)m_jpsiraw_lorentzvector->At(i))->Pz();
        double nrj = ((TLorentzVector*)m_jpsiraw_lorentzvector->At(i))->E();
        double m = sqrt( nrj*nrj - px*px - py*py - pz*pz );

        double interx = pow( (px/m)/sigmax, 2.);
        double intery = pow( (py/m)/sigmay, 2.);
        double interz = pow( (pz/m)/sigmaz, 2.);

        m_jpsiraw_sigmaL3D[i] = pow( interx + intery + interz , -0.5);
        //std::cout << "sigmaL3D = " << m_jpsiraw_sigmaL3D[i] << std::endl;

        double part1 = 0.;
        double part2 = 0.;
        double part3 = 0.;
        if (vtx.size() >0) {
          part1 = (px/m)*pow(m_jpsiraw_sigmaL3D[i]/sigmax,2.)*( svPos.x() - vtx[0].x());
          part2 = (py/m)*pow(m_jpsiraw_sigmaL3D[i]/sigmay,2.)*( svPos.y() - vtx[0].y());
          part3 = (pz/m)*pow(m_jpsiraw_sigmaL3D[i]/sigmaz,2.)*( svPos.z() - vtx[0].z());
        }  else {
          part1 = (px/m)*pow(m_jpsiraw_sigmaL3D[i]/sigmax,2.)*svPos.x();
          part2 = (py/m)*pow(m_jpsiraw_sigmaL3D[i]/sigmay,2.)*svPos.y();
          part3 = (pz/m)*pow(m_jpsiraw_sigmaL3D[i]/sigmaz,2.)*svPos.z();
        }

        m_jpsiraw_L3D[i] = fabs(part1 + part2 + part3);
        //std::cout << "L3D = " << m_jpsiraw_L3D[i] << std::endl;  

        m_jpsiraw_L3DoverSigmaL3D[i] = m_jpsiraw_L3D[i]/m_jpsiraw_sigmaL3D[i];
        //std::cout << "(L/sigma)3D = " << m_jpsiraw_L3DoverSigmaL3D[i] << std::endl;    

        if ( !vtxHandle.isValid() ) {
          std::cout << "PFpartExtractor::writeInfo(): vtxHandle is not valid..." << std::endl;

        }

        vector< RefCountedKinematicParticle > jpsi1_children = vertexFitTree->finalStateParticles();
        if ( jpsi1_children.size() != 2 ) {
          cout << " Warning Jpsi1 children size not equal to 2..." << endl;
        } else {
          // Order is : x,y,z,px,py,pz,m
          AlgebraicVector7 par1 = jpsi1_children[0]->currentState().kinematicParameters().vector();
          double e1 = jpsi1_children[0]->currentState().kinematicParameters().energy();
          new((*m_jpsiraw_mu1_lorentzvector)[i]) TLorentzVector(par1(3),par1(4),par1(5),e1);

          AlgebraicVector7 par2 = jpsi1_children[1]->currentState().kinematicParameters().vector();
          double e2 = jpsi1_children[1]->currentState().kinematicParameters().energy();
          new((*m_jpsiraw_mu2_lorentzvector)[i]) TLorentzVector(par2(3),par2(4),par2(5),e2);
        }

      } else cout << "Decay vertex Not valid\n";

    }
  }

  fillTree();
}


// Method initializing everything (to do for each event)

void PFpartExtractor::reset()
{

  // PF particle tree

  m_size = 0;
  m_pf_size = 0;

  for (int i=0;i<m_pfpart_MAX;++i) {
    m_pf_vx[i]         = 0.;
    m_pf_vy[i]         = 0.;
    m_pf_vz[i]         = 0.;
    m_pf_charge[i]     = 0;
    m_pf_pdgid[i]      = 0.;
    m_pf_trkLayer[i]   = 0;
    m_pf_pixLayer[i]   = 0;
    m_pf_trknormChi2[i]= 0.;
    m_pf_numberOfChambers[i]        = 0;
    m_pf_numberOfMatchedStations[i] = 0;
    m_pf_isGlobalMuon[i] = false;
    m_pf_isTrackerMuon[i] = false;
    m_pf_isStandAloneMuon[i] = false;
    m_pf_isCaloMuon[i] = false;
    m_pf_isPFMuon[i] = false;
    m_pf_isRPCMuon[i] = false;
  }

  m_pf_lorentzvector->Clear();

  // Jpsi tree
  m_jpsimu_lorentzvector->Clear();

  m_jpsiraw_lorentzvector->Clear();
  m_jpsiraw_mu1_lorentzvector->Clear();
  m_jpsiraw_mu2_lorentzvector->Clear();

  m_jpsi_size = 0;

  for (int i=0;i<m_jpsi_MAX;++i) {
    m_jpsi_indpf1[i] = 0;
    m_jpsi_indpf2[i] = 0;

    m_jpsiraw_vx[i] = 0.;
    m_jpsiraw_vy[i] = 0.;
    m_jpsiraw_vz[i] = 0.;
    m_jpsiraw_vtxvalid[i] = false;
    m_jpsiraw_vtxchi2[i] = 0.;
    m_jpsiraw_ndf[i] = 0.;

    m_jpsiraw_L3D[i] = 0.;
    m_jpsiraw_sigmaL3D[i] = 0.;
    m_jpsiraw_L3DoverSigmaL3D[i] = 0.;

  }

}


void PFpartExtractor::fillTree()
{
  m_tree_pfpart->Fill();
}

//
// Method getting the info from an input file
//

void PFpartExtractor::getInfo(int ievt) 
{
  m_tree_pfpart->GetEntry(ievt); 
}

