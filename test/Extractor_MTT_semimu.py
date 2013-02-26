from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()

options.register ('globalTag',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "The globaltag to use")

options.parseArguments()
if len(options.globalTag) == 0:
  raise Exception("You _must_ pass a globalTag options to this script. Use --help for more informations")

from Extractor_MTT_common import *

process = createExtractorProcess(False, True, useShiftCorrectedMET = True, globalTag = options.globalTag)

process.source.fileNames = cms.untracked.vstring(
    '/store/user/sbrochet/SingleMu/SingleMu_Run2012B-TOPMuPlusJets-13Jul2012_22Nov12-v1/482cfd8beeb5bd50ce95db8c4b04846e/patTuple_88_1_jP1.root'
    )

process.maxEvents.input = 5000
