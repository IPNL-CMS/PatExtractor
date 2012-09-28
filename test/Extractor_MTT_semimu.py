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
process.GlobalTag.globaltag = 'GR_R_53_V13::All'
process.PATextraction.doMC  = False

#Input PAT file to extract
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'/store/user/sbrochet/MuHad/MuHad_Run2012A_DCSONLY_v1/02979db18d11879dfe0836140c8c76cf/patTuple_76_2_hm9.root'
      '/store/user/sbrochet/SingleMu/SingleMu_Run2012B-TOPMuPlusJets-PromptSkim_08June/ff2bcae921e303fb6bcdd4793030d79d/patTuple_64_1_Nvx.root'      
      ),                           
    duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
    )

#Output extracted file name
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
process.PATextraction.jetCorrectorLabel = "ak5PFL1FastL2L3Residual"

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

    "W_mass                   80.399",
    "Top_mass                 172.0",
    "W_mass_err               10",
    "Top_mass_err             15.2",
    "b_mass                   4.67",

    "doSemiMu                 1",
    "doSyst                   0",
    "systvalue                1",
    "doUseBTaginChi2          1",
    "doChoiceWKF              0",

    "trigger                  ^HLT_IsoMu17_eta2p1_TriCentralPF(NoPU)?Jet[0-9]{0,3}(_[0-9]{0,3}){0,2}_v[0-9]{0,2}$"
    )


#########################################
#
# Launch the job
#
#########################################


process.p = cms.Path(process.PATextraction)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
