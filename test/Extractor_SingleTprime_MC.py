from SingleTprime_analysis_cfg import *

process = createExtractorProcess(True, False, useShiftCorrectedMET = False, globalTag = "START53_V27")

process.source.fileNames = cms.untracked.vstring(
   '/store/user/sperries/ZPrimeToTTJets_M1000GeV_W10GeV_TuneZ2star_8TeV-madgraph-tauola/Zprime_1000_Narrow_START53_V7A_04Dec12-v1/bd09b58f34b981e2c3ef3678b9b096ed/patTuple_2_1_P9c.root' 
    )
process.maxEvents.input = 5000
