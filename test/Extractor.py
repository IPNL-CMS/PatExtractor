#########################################
#
# Base macro for launching the PatExtractor
#
# The macro is for tests
#
#########################################

def readFile(file):
  return cms.untracked.string(open(file).read())

import FWCore.ParameterSet.Config as cms

process = cms.Process("PATextractor2")


#########################################
#
# Main configuration statements
#
#########################################

isMC = True

process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Extractors.PatExtractor.PAT_extractor_cff")

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1) #
  )

#Global tag and data type choice
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
process.PATextraction.isMC  = isMC
process.PATextraction.doMC  = isMC

#Input PAT file to extract
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
      'file:Zprime750_13TeV_miniaod.root'
      ),
  skipEvents=cms.untracked.uint32(15278),
  duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' ),
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
process.PATextraction.MET_PF.input  = cms.InputTag("patMETsPFlow")

process.PATextraction.doVertex   = True
process.PATextraction.vtx_tag    = cms.InputTag( "goodOfflinePrimaryVertices" )
process.PATextraction.doHLT      = True

if not isMC:
  process.PATextraction.triggersXML = readFile("triggers.xml")

# Jets correction : needs a valid global tags, or an external DB where JEC are stored
process.PATextraction.jet_PF.redoJetCorrection = False

if isMC:
  process.PATextraction.jet_PF.jetCorrectorLabel = "ak5PFchsL1FastL2L3"
else:
  process.PATextraction.jet_PF.jetCorrectorLabel = "ak5PFchsL1FastL2L3Residual"

process.PATextraction.jet_PF.doJER = False # Disable automatically on data

# JER systematics:
# Use -1 for 1-sigma down, 0 for nominal correction, and 1 for 1-sigma up
process.PATextraction.jet_PF.jerSign = 0

# JES systematics:
# Use -1 for 1-sigma down, 0 for nominal correction, and 1 for 1-sigma up
process.PATextraction.jet_PF.jesSign = 0
process.PATextraction.jet_PF.jes_uncertainties_file = cms.untracked.string("Extractors/PatExtractor/data/START53_V23_Uncertainty_AK5PFchs.txt")

process.PATextraction.MET_PF.redoMetPhiCorrection   = False
process.PATextraction.MET_PF.redoMetTypeICorrection = False # Automatically true if redoJetCorrection is True



#########################################
#
# Launch the job
#
#########################################

process.p = cms.Path(process.PATextraction)
process.MessageLogger.cerr.FwkReport.reportEvery = 1
