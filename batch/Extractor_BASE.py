#########################################
#
# Base macro for launching the PatExtractor
#
# The macro is sent via the script:
# /batch/runPAT_extraction.csh
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
process.load("Extractors.PatExtractor_2.PAT_extractor_cff")

#Other statements

process.GlobalTag.globaltag = 'GLOBALTAGNAME'

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) #
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    INPUTFILENAME
    ),                           
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

process.PATextraction.doMC = MCORNOTMC
process.PATextraction.extractedRootFile=cms.string('extracted.root')



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


process.PATextraction.doMuon     = True
process.PATextraction.doElectron = True
process.PATextraction.doJet      = True
process.PATextraction.doMET      = True
process.PATextraction.doVertex   = True


#########################################
#
# Launch the job
#
#########################################


process.p = cms.Path(process.PATextraction)
process.MessageLogger.cerr.FwkReport.reportEvery = 100
