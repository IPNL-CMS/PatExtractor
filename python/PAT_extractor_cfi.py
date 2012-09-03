import FWCore.ParameterSet.Config as cms

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
## Then you define the content of the output file (all set to false, turn on on request in you config files)
##
                               
   # Add HLT information
   doHLT         = cms.untracked.bool(False),
                               
   # Add MC information
   doMC          = cms.untracked.bool(False),
   MC_tag        = cms.InputTag( "" ),
                               
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
   jet_tag       = cms.InputTag( "selectedPatJetsPFlow" ),

   # Jets correction : needs a valid global tags, or an external DB where JEC are stored
   correctJets        = cms.untracked.bool(False),
   jetCorrectorLabel  = cms.untracked.string("ak5PFL1FastL2L3"), # Use "ak5PFL1FastL2L3" for MC and "ak5PFL1FastL2L3Residual" for Data

   # Add MET information
   doMET         = cms.untracked.bool(False),
   met_tag       = cms.InputTag( "patMETsPFlow" ),

   # Add PV information
   doVertex      = cms.untracked.bool(False),
   vtx_tag       = cms.InputTag( "offlinePrimaryVertices" ),

   # Add Track information
   doTrack       = cms.untracked.bool(False),
   trk_tag       = cms.InputTag( "generalTracks" ),


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
