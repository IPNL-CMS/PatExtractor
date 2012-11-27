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

process = createExtractorProcess(False, False, options.globalTag)

process.source.fileNames = cms.untracked.vstring(
    '/store/user/sbrochet/SingleElectron/SingleElectron_Run2012B-TOPElePlusJets-PromptSkim_part2_16June/ff2bcae921e303fb6bcdd4793030d79d/patTuple_44_1_U75.root'
    )
