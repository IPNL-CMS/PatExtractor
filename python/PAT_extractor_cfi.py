import os

import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.JetCorrectionProducers_cff import *

# Do some CHS stuff
ak5PFchsL1Fastjet  = ak5PFL1Fastjet.clone(algorithm = 'AK5PFchs')
ak5PFchsL2Relative = ak5PFL2Relative.clone(algorithm = 'AK5PFchs')
ak5PFchsL3Absolute = ak5PFL3Absolute.clone(algorithm = 'AK5PFchs')
ak5PFchsResidual   = ak5PFResidual.clone(algorithm = 'AK5PFchs')
ak5PFchsL1FastL2L3 = cms.ESProducer(
  'JetCorrectionESChain',
  correctors = cms.vstring('ak5PFchsL1Fastjet', 'ak5PFchsL2Relative','ak5PFchsL3Absolute')
)
ak5PFchsL1FastL2L3Residual = cms.ESProducer(
  'JetCorrectionESChain',
  correctors = cms.vstring('ak5PFchsL1Fastjet', 'ak5PFchsL2Relative','ak5PFchsL3Absolute', 'ak5PFchsResidual')
)

from Extractor_ScaleFactors import loadMuonScaleFactor, loadElectronScaleFactor, loadLightJetsScaleFactor
from RecoBTag.PerformanceDB.BTagPerformanceDBWinter13 import *

rootPath = os.path.join(os.environ["CMSSW_BASE"], "src/Extractors/PatExtractor/python")

PATextraction = cms.EDAnalyzer("PatExtractor",




##
## First you define the name of the  output ROOTfile
##
                               
  extractedRootFile = cms.string('extracted.root'),


##
## Then the name of the input ROOTfile, if you start from already extracted file
##
                               
  inputRootFile     = cms.string('default.root'),


##
## Here you tell is you start from a PATuple (True) or an extracted ROOTuple (False)
##

  fillTree = cms.untracked.bool(True),


##
## Are you running on data or MC?
##
   isMC          = cms.untracked.bool(True),


##
## Then you define the content of the output file (all set to false, turn on on request in you config files)
##
                               
   # Add HLT information
   doHLT         = cms.untracked.bool(False),
                               
   # Add MC information
   doMC          = cms.untracked.bool(False),
   MC_tag        = cms.InputTag( "" ),
   doMCjpsi      = cms.untracked.bool(False),
                             
   # Add Photon information
   doPhoton      = cms.untracked.bool(False),
   photon_tag    = cms.InputTag( "selectedPatPhotons" ),

   # Add Electron information
   doElectron    = cms.untracked.bool(False),
   electron_tag  = cms.InputTag( "selectedPatElectronsPFlow" ),

   # Add Muon information
   doMuon        = cms.untracked.bool(False),
   muon_tag      = cms.InputTag( "selectedPatMuonsPFlow" ),

   # Add Jet information
   doJet         = cms.untracked.bool(False),
   jet_PF        = cms.PSet(
       input              = cms.InputTag("selectedPatJetsPFlow"),

       # Jets correction : needs a valid global tags, or an external DB where JEC are stored
       redoJetCorrection      = cms.untracked.bool(False),
       jetCorrectorLabel      = cms.string("ak5PFchsL1FastL2L3"), # Use "ak5PFchsL1FastL2L3" for MC and "ak5PFchsL1FastL2L3Residual" for Data
       doJER                  = cms.untracked.bool(True),
       jerSign                = cms.untracked.int32(0), # Use 0 for no JER systematics, 1 for 1-sigma up and -1 for 1-sigma down
       jesSign                = cms.untracked.int32(0), # Use 0 for no JES systematics, 1 for 1-sigma up and -1 for 1-sigma down
       jes_uncertainties_file = cms.untracked.string(""),
       doLooseJetID           = cms.untracked.bool(True),
       useGlobalTagForJEC     = cms.untracked.bool(True),
       jecPayload             = cms.untracked.string("Extractors/PatExtractor/data/jec_payloads.xml"), 
       jecJetAlgo             = cms.untracked.string("AK5PFchs") 
   ),

   # Add MET information
   doMET         = cms.untracked.bool(False),
   MET_PF        = cms.PSet(
       input                  = cms.InputTag("patMETsPFlow"),

       redoMetPhiCorrection   = cms.untracked.bool(False),
       redoMetTypeICorrection = cms.untracked.bool(False),
       saveUnclusteredParticles = cms.untracked.bool(False)
   ),

   # Add KVF information (for J/psi and D0 reconstruction -- will use jet_PF infos)
   doKVF         = cms.untracked.bool(False), 
   jpsi_KVF      = cms.PSet(
        muJpsiMinPt   = cms.untracked.double(4.),
        jpsiMassMin   = cms.untracked.double(2.8),
        jpsiMassMax   = cms.untracked.double(3.4),
   ),
   d0_KVF      = cms.PSet(
        nTrD0Max      = cms.untracked.uint32(3),
        trD0MinPt     = cms.untracked.double(0.2),
   ),

   # Add PV information
   doVertex      = cms.untracked.bool(False),
   vtx_tag       = cms.InputTag( "offlinePrimaryVertices" ),

   # Add Track information
   doTrack       = cms.untracked.bool(False),
   trk_tag       = cms.InputTag( "generalTracks" ),

   # Add Track information
   doPF          = cms.untracked.bool(False),
   pf_tag        = cms.InputTag( "particleFlow" ),

   # Scale factors
   muon_scale_factors_tighteff_tightiso = loadMuonScaleFactor(os.path.join(rootPath, "MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl"), os.path.join(rootPath, "MuonEfficiencies_Run2012ReReco_53X.pkl"), "Tight", "combRelIsoPF04dBeta<012_Tight"),
   muon_scale_factors_tighteff_looseiso = loadMuonScaleFactor(os.path.join(rootPath, "MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl"), os.path.join(rootPath, "MuonEfficiencies_Run2012ReReco_53X.pkl"), "Tight", "combRelIsoPF04dBeta<02_Tight"),
   muon_scale_factors_looseeff_looseiso = loadMuonScaleFactor(os.path.join(rootPath, "MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl"), os.path.join(rootPath, "MuonEfficiencies_Run2012ReReco_53X.pkl"), "Loose", "combRelIsoPF04dBeta<02_Loose"),
   muon_scale_factors = cms.vstring("muon_scale_factors_looseeff_looseiso", "muon_scale_factors_tighteff_looseiso", "muon_scale_factors_tighteff_tightiso"),

   electron_scale_factors_tighteff_tightiso = loadElectronScaleFactor(os.path.join(rootPath, "Electron_scale_factors.json"), os.path.join(rootPath, "Electrons_ScaleFactors_Reco_8TeV.root"), "tight"),
   electron_scale_factors_looseeff_tightiso = loadElectronScaleFactor(os.path.join(rootPath, "Electron_scale_factors.json"), os.path.join(rootPath, "Electrons_ScaleFactors_Reco_8TeV.root"), "loose"),
   electron_scale_factors = cms.vstring("electron_scale_factors_tighteff_tightiso", "electron_scale_factors_looseeff_tightiso"),

   b_tagging_scale_factors_b_jets = cms.PSet(
      jet_type = cms.string("b"),
      from_globaltag = cms.bool(True),
      payload = cms.string("MUJETSWPBTAGNOTTBARCSVM")
      ),
   b_tagging_scale_factors_c_jets = cms.PSet(
      jet_type = cms.string("c"),
      from_globaltag = cms.bool(True),
      payload = cms.string("MUJETSWPBTAGNOTTBARCSVM")
      ),
   b_tagging_scale_factors_light_jets = cms.PSet(
      jet_type = cms.string("light"),
      from_globaltag = cms.bool(False),
      scale_factors = loadLightJetsScaleFactor()
      ),
   b_tagging_scale_factors = cms.vstring("b_tagging_scale_factors_b_jets", "b_tagging_scale_factors_c_jets", "b_tagging_scale_factors_light_jets"),


##
## Finally you put some details on the analysis
##

   doDimuon      = cms.untracked.bool(False),
   do4TopHLT     = cms.untracked.bool(False),
   n_events = cms.untracked.int32(10000),  # How many events you want to analyze (only if fillTree=False)

   # The analysis settings could be whatever you want
   # 
   # Format is "STRING VALUE" where STRING is the name of the cut, and VALUE the value of the cut

   # Here we define for example the cuts for the dimuon analysis
                               
   analysisSettings = cms.untracked.vstring()
)
