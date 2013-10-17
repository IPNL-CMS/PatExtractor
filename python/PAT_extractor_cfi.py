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
   isMC          = cms.untracked.bool(False),


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
       redoJetCorrection  = cms.untracked.bool(False),
       jetCorrectorLabel  = cms.string("ak5PFchsL1FastL2L3"), # Use "ak5PFchsL1FastL2L3" for MC and "ak5PFchsL1FastL2L3Residual" for Data
       doJER              = cms.untracked.bool(True),
       jerSign            = cms.untracked.int32(0), # Use 0 for no JER systematics, 1 for 1-sigma up and -1 for 1-sigma down
       jesSign            = cms.untracked.int32(0), # Use 0 for no JES systematics, 1 for 1-sigma up and -1 for 1-sigma down
   ),

   # Add MET information
   doMET         = cms.untracked.bool(False),
   MET_PF        = cms.PSet(
       input                  = cms.InputTag("patMETsPFlow"),

       redoMetPhiCorrection   = cms.untracked.bool(False),
       redoMetTypeICorrection = cms.untracked.bool(False),
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

##
## Finally you put some details on the analysis
##

   doMtt         = cms.untracked.bool(False),
   doDimuon      = cms.untracked.bool(False),
   do4TopHLT     = cms.untracked.bool(False),
   n_events = cms.untracked.int32(10000),  # How many events you want to analyze (only if fillTree=False)

   # The analysis settings could be whatever you want
   # 
   # Format is "STRING VALUE" where STRING is the name of the cut, and VALUE the value of the cut

   # Here we define for example the cuts for the dimuon analysis
                               
   analysisSettings = cms.untracked.vstring()
)
