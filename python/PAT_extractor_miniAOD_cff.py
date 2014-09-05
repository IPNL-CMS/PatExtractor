import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.MagneticField_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducers_cff import *

from Extractors.PatExtractor.PAT_extractor_cfi import *

# Customize input collection to use ones from miniAOD
PATextraction.extractors.MC.parameters.input = cms.InputTag("prunedGenParticles")

PATextraction.extractors.HLT.parameters.prescales = cms.InputTag("patTrigger")

PATextraction.extractors.track.enable = False

PATextraction.extractors.PFpart.parameters.input = cms.InputTag('packedPFCandidates')
PATextraction.extractors.PFpart.parameters.vertices = cms.InputTag('offlineSlimmedPrimaryVertices')

PATextraction.extractors.Vertices.parameters.input = cms.InputTag('offlineSlimmedPrimaryVertices')

PATextraction.extractors.electron_PF.parameters.input = cms.InputTag("slimmedElectrons")
PATextraction.extractors.electron_PF.parameters.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
PATextraction.extractors.electron_PF.parameters.conversions = cms.InputTag("reducedEgamma", "reducedConversions", "PAT")
PATextraction.extractors.electron_PF.parameters.rho = cms.InputTag("fixedGridRhoAll")

PATextraction.extractors.muon_PF.parameters.input = cms.InputTag("slimmedMuons")
PATextraction.extractors.muon_PF.parameters.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")

PATextraction.extractors.jetmet.parameters.input_jets = cms.InputTag("slimmedJets")
PATextraction.extractors.jetmet.parameters.input_met = cms.InputTag("slimmedMETs")
PATextraction.extractors.jetmet.parameters.input_raw_met = cms.InputTag("NotAvailable")
PATextraction.extractors.jetmet.parameters.pf_candidates = cms.InputTag("packedPFCandidates")
PATextraction.extractors.jetmet.parameters.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
PATextraction.extractors.jetmet.parameters.rho = cms.InputTag("fixedGridRhoAll")


PATextraction.extractors.photon.enable = False
PATextraction.extractors.photon.parameters.input = cms.InputTag("slimmedPhotons")
