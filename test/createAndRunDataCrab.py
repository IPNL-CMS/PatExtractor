#! /usr/bin/env python

import os, datetime, pwd

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
(options, args) = parser.parse_args()

datasets = [

    ["/ElectronHad/sbrochet-ElectronHad_Run2012A-13Jul2012_22Nov12-v1-482cfd8beeb5bd50ce95db8c4b04846e/USER", "ElectronHad_Run2012A-13Jul2012", "FT_53_V6_AN2"],
    ["/ElectronHad/sbrochet-ElectronHad_Run2012A-recover-06Aug2012_22Nov12-v1-a4e328caacc917c09c291e66055c015f/USER", "ElectronHad_Run2012A-recover-06Aug2012", "FT_53_V6_AN2"],
    ["/SingleElectron/sbrochet-SingleElectron_Run2012B-TOPElePlusJets-13Jul2012_22Nov12-v1-482cfd8beeb5bd50ce95db8c4b04846e/USER", "SingleElectron_Run2012B-TOPElePlusJets-13Jul2012", "FT_53_V6_AN2"],
    ["/SingleElectron/sbrochet-SingleElectron_Run2012C-TOPElePlusJets-24Aug2012_22Nov12-v1-2648d2ac41c8a47bf1ce16deec221c74/USER", "SingleElectron_Run2012C-TOPElePlusJets-24Aug2012", "FT_53_V10_AN2"],
    ["/SingleElectron/sbrochet-SingleElectron_Run2012C-TOPElePlusJets-PromptSkim_22Nov12-v1-a5e8ad198ec0ef3c5777633fcf11cc8c/USER", "SingleElectron_Run2012C-TOPElePlusJets-PromptSkim", "GR_P_V41_AN2"],
    ["/SingleElectron/sbrochet-SingleElectron_Run2012D-TOPElePlusJets-PromptSkim_22Nov12-v1-ba18d001f77ddfaed5a189360146b128/USER", "SingleElectron_Run2012D-TOPElePlusJets-PromptSkim", "GR_P_V42_AN2"],

    ["/MuHad/sbrochet-MuHad_Run2012A-13Jul2012_22Nov12-v1-482cfd8beeb5bd50ce95db8c4b04846e/USER", "MuHad_Run2012A-13Jul2012", "FT_53_V6_AN2"],
    ["/MuHad/sbrochet-MuHad_Run2012A-recover-06Aug2012_22Nov12-v1-a4e328caacc917c09c291e66055c015f/USER", "MuHad_Run2012A-recover-06Aug2012", "FT_53_V6C_AN2"],
    ["/SingleMu/sbrochet-SingleMu_Run2012B-TOPMuPlusJets-13Jul2012_22Nov12-v1-482cfd8beeb5bd50ce95db8c4b04846e/USER", "SingleMu_Run2012B-TOPMuPlusJets-13Jul2012", "FT_53_V6_AN2"],
    ["/SingleMu/sbrochet-SingleMu_Run2012C-TOPMuPlusJets-24Aug2012_22Nov12-v1-2648d2ac41c8a47bf1ce16deec221c74/USER", "SingleMu_Run2012C-TOPMuPlusJets-24Aug2012", "FT_53_V10_AN2"],
    ["/SingleMu/sbrochet-SingleMu_Run2012C-TOPMuPlusJets-PromptSkim_22Nov12-v1-a5e8ad198ec0ef3c5777633fcf11cc8c/USER", "SingleMu_Run2012C-TOPMuPlusJets-PromptSkim", "GR_P_V41_AN2"],
    ["/SingleMu/sbrochet-SingleMu_Run2012D-TOPMuPlusJets-PromptSkim_22Nov12-v1-ba18d001f77ddfaed5a189360146b128/USER", "SingleMu_Run2012D-TOPMuPlusJets-PromptSkim", "GR_P_V42_AN2"]

    ]

# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

d = datetime.datetime.now().strftime("%d%b")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

for dataset in datasets:
  dataset_path = dataset[0]
  dataset_name = dataset[1]
  dataset_globaltag = dataset[2]

  ui_working_dir = ("crab_data_%s_%s") % (dataset_name, d)
  output_file = "crab_data_%s_%s.cfg" % (dataset_name, d)
  output_dir = ("MTT/Extracted/data/%s/%s" % (d, dataset_name))

  python_config = "";
  if "Electron" in dataset_path:
    python_config = "Extractor_MTT_semie.py"
  else:
    python_config = "Extractor_MTT_semimu.py"

  print("Creating config file for '%s'" % (dataset_path))
  print("\tName: %s" % dataset_name)
  print("\tGlobal tag: %s" % dataset_globaltag)
  print("\tOutput directory: %s" % output_dir)
  print("")
    

  os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@globaltag@#%s#g\" -e \"s#@pset@#%s#g\" crab_data.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir, email, dataset_globaltag, python_config, output_file))

  if options.run:
    cmd = "crab -create -submit -cfg %s" % (output_file)
    os.system(cmd);
