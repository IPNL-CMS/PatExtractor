from Extractor_MTT_common_jpsi import *

process = createExtractorProcess(False, True, useShiftCorrectedMET = True, globalTag = "FT_53_V21_AN6")

# To build Transient Tracks
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.source.fileNames = cms.untracked.vstring(
    'file:/gridgroup/cms/bouvier/data/test_extractor/patTuple_1_1_Uqw.root'
    )
process.maxEvents.input = -1
