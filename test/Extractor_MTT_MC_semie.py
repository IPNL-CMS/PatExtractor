from Extractor_MTT_common import *

process = createExtractorProcess(True, False, "START53_V7F")

process.source.fileNames = cms.untracked.vstring(
    '/store/user/sperries/ZPrimeToTTJets_M750GeV_W7p5GeV_TuneZ2star_8TeV-madgraph-tauola/Zprime_750_Narrow_2012_PF2PAT_v1/165778d6ec003db3c40b0ea37fd1f4fc/patTuple_1_1_pTr.root'
    )
