from Extractor_SingleTprime_common import *

process = createExtractorProcess(True, False, useShiftCorrectedMET = False, globalTag = "START53_V15A")

# To build Transient Tracks
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#process.source.fileNames = cms.untracked.vstring(
    #'file:/gridgroup/cms/bouvier/data/test_extractor/patTuple_1_1_Uqw.root'
    #'file:/gridgroup/cms/bouvier/data/test_extractor/patTuple_10_1_Wtz.root'
#    'file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/PatExtractor/Small_Sample_patTuple.root'
    #'file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/PatExtractor/test/extracted_mc.root'
#    )

process.maxEvents.input = -1
