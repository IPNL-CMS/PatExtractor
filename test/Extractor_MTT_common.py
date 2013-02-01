#########################################
#
# Base macro for launching the PatExtractor
#
# The macro is for tests
#
#########################################


import FWCore.ParameterSet.Config as cms

def readFile(file):
  return cms.untracked.string(open(file).read())

def createExtractorProcess(isMC, isSemiMu, useShiftCorrectedMET, globalTag):
  process = cms.Process("PATextractor2")


  #########################################
  #
  # Main configuration statements
  #
  #########################################

  process.load('Configuration/StandardSequences/Services_cff')
  process.load('Configuration/StandardSequences/GeometryIdeal_cff')
  process.load('Configuration/StandardSequences/MagneticField_38T_cff')
  process.load('Configuration/StandardSequences/EndOfProcess_cff')
  process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
  process.load("FWCore.MessageLogger.MessageLogger_cfi")
  process.load("Extractors.PatExtractor.PAT_extractor_cff")

  process.maxEvents = cms.untracked.PSet(
      input = cms.untracked.int32(10) #
      )

  #Global tag and data type choice
  process.GlobalTag.globaltag = '%s::All' % globalTag
  process.PATextraction.isMC  = isMC
  process.PATextraction.doMC  = isMC

  #Input PAT file to extract
  process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(
        ),                           
      duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
      )

  #Output extracted file name
  if isMC:
    process.PATextraction.extractedRootFile = cms.string('extracted_mc.root')
  else:
    process.PATextraction.extractedRootFile = cms.string('extracted.root')

  #########################################
  #
  # PAT extractor main options statements
  #
  #########################################

  #
  # Adapt it to your needs
  #
  # If you are lost, see the example here (PART 3.2):
  # http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.PHYTuto
  #
  # Here we just extract, and don't perform any analysis

  process.PATextraction.doMuon     = True
  process.PATextraction.doElectron = True
  process.PATextraction.doJet      = True

  process.PATextraction.doMET      = True
  if useShiftCorrectedMET:
    process.PATextraction.met_tag  = cms.InputTag("patMETsShiftCorrPFlow")
  else:
    process.PATextraction.met_tag  = cms.InputTag("patMETsPFlow")

  process.PATextraction.doVertex   = True
  process.PATextraction.vtx_tag    = cms.InputTag( "goodOfflinePrimaryVertices" )
  process.PATextraction.doHLT      = True

  if not isMC:
    if isSemiMu:
      process.PATextraction.triggersXML = readFile("triggers_mu.xml")
    else:
      process.PATextraction.triggersXML = readFile("triggers_e.xml")

  process.PATextraction.doMtt      = True

  # Jets correction : needs a valid global tags, or an external DB where JEC are stored
  process.PATextraction.correctJets       = True

  process.PATextraction.correctSysShiftMet       = True

  if isMC:
    process.PATextraction.jetCorrectorLabel = "ak5PFchsL1FastL2L3"
  else:
    process.PATextraction.jetCorrectorLabel = "ak5PFchsL1FastL2L3Residual"

  process.PATextraction.redoTypeIMET      = False

  # Analysis cuts
  import sys
  sys.path.append('.')
  
  if isSemiMu:
    import Extractor_MTT_analysis_cuts_semimu as anaSettings
  else:
    import Extractor_MTT_analysis_cuts_semie as anaSettings
    
  process.PATextraction.analysisSettings = anaSettings.analysisSettings

  # MTT analysis configuration
  process.PATextraction.mtt = cms.PSet(
      # ------------------------------------------------
      # settings for the KinFitter
      # ------------------------------------------------    
      maxNrIter = cms.uint32(500),
      maxDeltaS = cms.double(5e-05),
      maxF      = cms.double(0.0001),
      # ------------------------------------------------
      # select parametrisation
      # 0: EMom, 1: EtEtaPhi, 2: EtThetaPhi
      # ------------------------------------------------
      jetParametrisation = cms.uint32(1),
      lepParametrisation = cms.uint32(1),
      metParametrisation = cms.uint32(1),

      # ------------------------------------------------
      # set constraints
      # 1: Whad-mass, 2: Wlep-mass, 3: thad-mass,
      # 4: tlep-mass, 5: nu-mass, 6: equal t-masses
      # 7: sum-pt conservation
      # ------------------------------------------------
      constraints = cms.vuint32(1, 2, 3, 4),

      # ------------------------------------------------
      # set mass values used in the constraints
      # ------------------------------------------------
      mW   = cms.double(80.4),
      mTop = cms.double(173.),

      # ------------------------------------------------
      # set correction factor(s) for the jet energy resolution:
      # - (optional) eta dependence assumed to be symmetric
      #   around eta=0, i.e. parametrized in |eta|
      # - any negative value as last bin edge is read as "inf"
      # - make sure that number of entries in vector with
      #   bin edges = number of scale factors + 1
      # ------------------------------------------------
      jetEnergyResolutionScaleFactors = cms.vdouble(1.0),
      jetEnergyResolutionEtaBinning = cms.vdouble(0.0,-1.0))

  #########################################
  #
  # Launch the job
  #
  #########################################


  process.p = cms.Path(process.PATextraction)
  process.MessageLogger.cerr.FwkReport.reportEvery = 1000

  return process
