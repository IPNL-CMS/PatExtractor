import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.MagneticField_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

from Extractors.PatExtractor.PAT_extractor_cfi import *

#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
## turn on VID producer, indicate data format  to be
## DataFormat.AOD or DataFormat.MiniAOD, as appropriate

## Set up input/output depending on the format
## You can list here either AOD or miniAOD files, but not both types mixed
##
#useAOD = False

#if useAOD == True :
    #dataFormat = DataFormat.AOD
#else :
    #dataFormat = DataFormat.MiniAOD
 
#switchOnVIDPhotonIdProducer(process, dataFormat)
 
## define which IDs we want to produce
#my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']
 
##add them to the VID producer
#for idmod in my_id_modules:
    #setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

# Customize input collection to use ones from miniAOD
PATextraction.extractors.MC.parameters.input = cms.InputTag("prunedGenParticlesPFlow")

PATextraction.extractors.HLT.parameters.prescales = cms.InputTag("patTrigger")

PATextraction.extractors.track.enable = False

PATextraction.extractors.PFpart.parameters.input = cms.InputTag('packedPFCandidatesPFlow')
PATextraction.extractors.PFpart.parameters.vertices = cms.InputTag('offlineSlimmedPrimaryVertices')

PATextraction.extractors.Vertices.parameters.input = cms.InputTag('offlineSlimmedPrimaryVertices')

PATextraction.extractors.electron_PF.parameters.input = cms.InputTag("slimmedElectronsPFlow")
PATextraction.extractors.electron_PF.parameters.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
PATextraction.extractors.electron_PF.parameters.conversions = cms.InputTag("reducedEgamma", "reducedConversions", "PAT")
PATextraction.extractors.electron_PF.parameters.rho = cms.InputTag("fixedGridRhoAll")

PATextraction.extractors.electrons_loose.enable = cms.bool(True)
PATextraction.extractors.electrons_loose.parameters.input = cms.InputTag("slimmedAllElectronsPFlow")
PATextraction.extractors.electrons_loose.parameters.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
PATextraction.extractors.electrons_loose.parameters.conversions = cms.InputTag("reducedEgamma", "reducedConversions", "PAT")
PATextraction.extractors.electrons_loose.parameters.rho = cms.InputTag("fixedGridRhoAll")

PATextraction.extractors.muon_PF.parameters.input = cms.InputTag("slimmedMuonsPFlow")
PATextraction.extractors.muon_PF.parameters.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")

PATextraction.extractors.muons_loose.enable = cms.bool(True)
PATextraction.extractors.muons_loose.parameters.input = cms.InputTag("slimmedAllMuonsPFlow")
PATextraction.extractors.muons_loose.parameters.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")

PATextraction.extractors.jetmet.parameters.input_jets = cms.InputTag("slimmedJetsPFlow")
PATextraction.extractors.jetmet.parameters.input_met = cms.InputTag("patMETsPFlow")
PATextraction.extractors.jetmet.parameters.input_raw_met = cms.InputTag("NotAvailable")
PATextraction.extractors.jetmet.parameters.pf_candidates = cms.InputTag("packedPFCandidatesPFlow")
PATextraction.extractors.jetmet.parameters.vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
PATextraction.extractors.jetmet.parameters.rho = cms.InputTag("fixedGridRhoAll")


PATextraction.extractors.photon.enable = False
PATextraction.extractors.photon.parameters.input = cms.InputTag("slimmedPhotonsPFlow")

#ExtractorSequence = cms.Sequence(egmPhotonsIDs + PATextraction)
