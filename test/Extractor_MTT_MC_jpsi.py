from Extractor_MTT_common_jpsi import *

process = createExtractorProcess(True, True, useShiftCorrectedMET = True, globalTag = "START53_V19PR")

# To build Transient Tracks
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.source.fileNames = cms.untracked.vstring(
    #'file:/gridgroup/cms/bouvier/data/test_extractor/patTuple_1_1_Uqw.root'
    'file:/gridgroup/cms/bouvier/data/test_extractor/patTuple_10_1_Wtz.root'
    )
process.maxEvents.input = -1
