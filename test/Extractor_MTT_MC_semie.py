#########################################
#
# Base macro for launching the PatExtractor
#
# The macro is for tests
#
#########################################


import FWCore.ParameterSet.Config as cms

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

process.options = cms.untracked.PSet(
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000) #
    )

#Global tag and data type choice
process.GlobalTag.globaltag = 'START53_V7A::All'
process.PATextraction.doMC  = True

#Input PAT file to extract
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/user/sbrochet/TTJets_TuneZ2star_8TeV-madgraph-tauola/TTJets_2012_v1/265c9c69c37a8e555f9b98fa1aae946f/patTuple_36_1_VWJ.root'
      ),                           
    duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
    )

#Output extracted file name
process.PATextraction.extractedRootFile=cms.string('extracted_mc.root')



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
process.PATextraction.doVertex   = True
process.PATextraction.doHLT      = True

process.PATextraction.doMtt      = True

# Jets correction : needs a valid global tags, or an external DB where JEC are stored
process.PATextraction.correctJets       = True
process.PATextraction.jetCorrectorLabel = "ak5PFchsL1FastL2L3"

# Analysis cuts
from Extractor_MTT_analysis_cuts_semie import *
process.PATextraction.analysisSettings = analysisSettings

#########################################
#
# Launch the job
#
#########################################


process.p = cms.Path(process.PATextraction)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
