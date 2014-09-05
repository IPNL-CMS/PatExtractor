import os

import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.JetCorrectionProducers_cff import *

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

        n_events = cms.untracked.int32(10000),  # How many events you want to analyze (only if fillTree=False)


        ##
        ## Are you running on data or MC?
        ##
        isMC          = cms.untracked.bool(True),


        ##
        ## Then you define the content of the output file (all set to false, turn on on request in you config files)
        ##

        extractors    = cms.PSet(

            event = cms.PSet(
                type = cms.string("event_extractor"),
                enable = cms.bool(True),
                parameters = cms.PSet(
                    pileup_summary = cms.InputTag("addPileupInfo"),
                    generator = cms.InputTag("generator")
                    )
                ),

            MC = cms.PSet(
                type = cms.string("mc_extractor"),
                enable = cms.bool(False),
                parameters = cms.PSet(
                    input = cms.InputTag("genParticles"),
                    do_jpsi = cms.bool(True)
                    )
                ),

            HLT = cms.PSet(
                type = cms.string("hlt_extractor"),
                enable = cms.bool(True),
                parameters = cms.PSet(
                    input = cms.InputTag("TriggerResults", "", "HLT"),
                    triggers = cms.untracked.string("")
                    )
                ),

            track = cms.PSet(
                type = cms.string("track_extractor"),
                enable = cms.bool(False),
                parameters = cms.PSet(
                    input = cms.InputTag("generalTracks")
                    )
                ),

            PFpart = cms.PSet(
                type = cms.string("pfparticle_extractor"),
                enable = cms.bool(False),
                parameters = cms.PSet(
                    input = cms.InputTag("particleFlow"),
                    vertices = cms.InputTag("offlinePrimaryVertices")
                    )
                ),

            Vertices = cms.PSet(
                type = cms.string("vertex_extractor"),
                enable = cms.bool(True),
                parameters = cms.PSet(
                    input = cms.InputTag("offlinePrimaryVertices")
                    )
                ),

            electron_PF = cms.PSet(
                type = cms.string("electron_extractor"),
                enable = cms.bool(True),
                parameters = cms.PSet(
                    input = cms.InputTag("selectedPatElectronsPFlow"),
                    vertices = cms.InputTag("offlinePrimaryVertices"),
                    conversions = cms.InputTag("allConversions"),
                    beamspot = cms.InputTag("offlineBeamSpot"),
                    rho = cms.InputTag("kt6PFJets", "rho", "RECO")
                    )
                ),

            electrons_loose = cms.PSet(
                    type = cms.string("electron_extractor"),
                    enable = cms.bool(False),
                    parameters = cms.PSet(
                        input = cms.InputTag("selectedPatElectronsLoosePFlow"),
                        vertices = cms.InputTag("offlinePrimaryVertices"),
                        conversions = cms.InputTag("allConversions"),
                        beamspot = cms.InputTag("offlineBeamSpot"),
                        rho = cms.InputTag("kt6PFJets", "rho", "RECO")
                        )
                    ),

            muon_PF = cms.PSet(
                    type = cms.string("muon_extractor"),
                    enable = cms.bool(True),
                    parameters = cms.PSet(
                        input = cms.InputTag("selectedPatMuonsPFlow"),
                        vertices = cms.InputTag("offlinePrimaryVertices")
                        )
                    ),

            muons_loose = cms.PSet(
                    type = cms.string("muon_extractor"),
                    enable = cms.bool(False),
                    parameters = cms.PSet(
                        input = cms.InputTag("selectedPatMuonsLoosePFlow"),
                        vertices = cms.InputTag("offlinePrimaryVertices")
                        )
                    ),

            jetmet = cms.PSet(
                    type = cms.string("jet_met_extractor"),
                    enable = cms.bool(True),
                    parameters = cms.PSet(
                        input_jets = cms.InputTag("selectedPatJetsPFlow"),
                        input_met = cms.InputTag("patMETsPFlow"),
                        input_raw_met = cms.InputTag("patPFMetPFlow"),
                        pf_candidates = cms.InputTag("particleFlow"),
                        vertices = cms.InputTag("offlinePrimaryVertices"),
                        rho = cms.InputTag("kt6PFJets", "rho", "RECO"),

                        tree_name_jets = cms.string("jet_PF"),
                        tree_name_met = cms.string("MET_PF"),

                        # Jets correction : needs a valid global tags, or an external DB where JEC are stored
                        redoJetCorrection      = cms.untracked.bool(False),
                        jetCorrectorLabel      = cms.string("ak4PFchsL1FastL2L3"), # Use "ak5PFchsL1FastL2L3" for MC and "ak5PFchsL1FastL2L3Residual" for Data
                        doJER                  = cms.untracked.bool(True),
                        jerSign                = cms.untracked.int32(0), # Use 0 for no JER systematics, 1 for 1-sigma up and -1 for 1-sigma down
                        jesSign                = cms.untracked.int32(0), # Use 0 for no JES systematics, 1 for 1-sigma up and -1 for 1-sigma down
                        jes_uncertainties_file = cms.untracked.string(""),
                        doLooseJetID           = cms.untracked.bool(True),
                        useGlobalTagForJEC     = cms.untracked.bool(True),
                        jecPayload             = cms.untracked.string("Extractors/PatExtractor/data/jec_payloads.xml"), 
                        jecJetAlgo             = cms.untracked.string("AK4PFchs"), 

                        redoMetPhiCorrection   = cms.untracked.bool(False),
                        redoMetTypeICorrection = cms.untracked.bool(False),
                        saveUnclusteredParticles = cms.untracked.bool(False)
                        )
                    ),

            photon = cms.PSet(
                    type = cms.string("photon_extractor"),
                    enable = cms.bool(False),
                    parameters = cms.PSet(
                        input = cms.InputTag("selectedPatPhotons"),
                        matched_electron = cms.InputTag("photonPFIsolation", "hasMatchedPromptElectron", "PAT"),
                        charged_hadrons_iso = cms.InputTag("photonPFIsolation", "chargedHadronsIsolation", "PAT"),
                        neutral_hadrons_iso = cms.InputTag("photonPFIsolation", "neutralHadronsIsolation", "PAT"),
                        photons_iso = cms.InputTag("photonPFIsolation", "photonsIsolation", "PAT"),
                        rho = cms.InputTag("kt6PFJets", "rho", "RECO")
                        )
                    )
            ),

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
)
