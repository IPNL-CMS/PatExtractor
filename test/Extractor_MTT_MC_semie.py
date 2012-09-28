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
process.GlobalTag.globaltag = 'START53_V11::All'
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
process.PATextraction.jetCorrectorLabel = "ak5PFL1FastL2L3"

# Analysis cuts
process.PATextraction.analysisSettings = cms.untracked.vstring(
    "VTX_Ndof_Min             4",
    
    "MET_Pt_Min               20",
    
    "MU_Pt_min_loose          10",
    "MU_Eta_max_loose         2.5",
    "MU_Iso_min               0.125",
    "MU_Pt_min                20",
    "MU_Eta_max               2.4",
    "MU_normChi2_max          10",
    "MU_nValTrackHits_min     10",
    "MU_nMatches_min          1",
    "MU_nValPixHits_min       1",
    "MU_dB_min                0.02",
    "MU_ePt_min               15",
    "MU_eEta_max              2.5",
    "MU_eEtaW_min             1.4442",
    "MU_eEtaW_max             1.5560",
    "MU_eIso_min              0.2",

    "ELE_Iso_min              0.1",
    "ELE_Pt_min               30",
    "ELE_Eta_max              2.5",
    "ELE_Zmass                91",
    "ELE_Zwin                 15",
    "ELE_dB_min               0.02",
    
    "JET_Pt_min               30",
    "JET_Eta_max              2.4",

    "JET_btag_CSVL_min        0.244",
    "JET_btag_CSVM_min        0.679",
    "JET_btag_CSVT_min        0.898",
    "JET_btag_TCHPT_min       3.41",
    "JET_btag_SSVHEM_min      1.74",
    "JET_btag_SSVHPT_min      2",

    "W_mass                   80.399",
    "Top_mass                 172.0",
    "W_mass_err               10",
    "Top_mass_err             15.2",
    "b_mass                   4.67",

    "doSemiMu                 0",
    "doSyst                   0",
    "systvalue                1",
    "doUseBTaginChi2          1",
    "doChoiceWKF              0",

    "trigger                  ^HLT_Ele25_CaloIdVT_CaloIso(VL|T)_TrkId(VL|T)_TrkIsoT_TriCentralPF(NoPU)?Jet30_v[0-9]{0,2}$"
    )


#########################################
#
# Launch the job
#
#########################################


process.p = cms.Path(process.PATextraction)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
