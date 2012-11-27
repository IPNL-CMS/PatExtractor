#! /usr/bin/env python

import os, copy, datetime, pwd, re

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
(options, args) = parser.parse_args()

datasets = {
    # Single anti-top
    "/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/sbrochet-Tbar_t-channel_2012_PF2PAT_v1-265c9c69c37a8e555f9b98fa1aae946f/USER": "Tbar_t-channel",
    "/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/sbrochet-Tbar_tW-channel_2012_PF2PAT_v1-265c9c69c37a8e555f9b98fa1aae946f/USER": "Tbar_tW-channel",
    "/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/sbrochet-Tbar_s-channel_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "Tbar_s-channel",
    
    # Single top
    "/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/sbrochet-T_tW-channel_2012_PF2PAT_v1-265c9c69c37a8e555f9b98fa1aae946f/USER": "T_tW-channel",
    "/T_s-channel_TuneZ2star_8TeV-powheg-tauola/sbrochet-T_s-channel_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "T_s-channel",

    ## TT + jets
    "/TTJets_TuneZ2star_8TeV-madgraph-tauola/sbrochet-TTJets_2012_v1-265c9c69c37a8e555f9b98fa1aae946f/USER": "TTJets_MassiveBinDECAY",
    
    "/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/sbrochet-DYJetsToLL_M-50_2012_PF2PAT_v1-265c9c69c37a8e555f9b98fa1aae946f/USER": "DYJetsToLL_M-50",
    "/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/sbrochet-WJetsToLNu_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "WJetsToLNu"

    #"/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_20_30_EMEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "QCD_Pt_20_30_EMEnriched",
    #"/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_30_80_EMEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "QCD_Pt_30_80_EMEnriched",
    #"/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_80_170_EMEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "QCD_Pt_80_170_EMEnriched",
    #"/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_20_MuEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "QCD_Pt_20_MuEnriched"
    }

# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

d = datetime.datetime.now().strftime("%d%b%y")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

from string import Template
multicrab = Template(r"""[MULTICRAB]
cfg=crab_MC.cfg.template.ipnl

[COMMON]
USER.ui_working_dir = ${ui_working_dir}
USER.eMail = ${email}
CMSSW.datasetpath = ${dataset}
"""
)

multicrab_semie = r"""
[${name}_semie]
CMSSW.pset = Extractor_MTT_MC_semie.py
USER.user_remote_dir = ${remote_dir_semie}
"""

multicrab_semimu = r"""
[${name}_semimu]
CMSSW.pset = Extractor_MTT_MC_semimu.py
USER.user_remote_dir = ${remote_dir_semimu}
"""

for dataset_path, dataset_name in datasets.items():

  #dataset_globaltag = re.search('START\d{0,2}_V\d[A-Z]?', dataset_path).group(0)

  #publish_name = "%s_%s_%s-v%d" % (dataset_name, dataset_globaltag, d, version)
  ui_working_dir = ("multicrab_MC_%s_%s") % (dataset_name, d)
  output_file = "multicrab_MC_%s_%s.cfg" % (dataset_name, d)

  output_dir_semie = ("MTT/Extracted/MC/Summer12/%s/semie/%s" % (d, dataset_name))
  output_dir_semimu = ("MTT/Extracted/MC/Summer12/%s/semimu/%s" % (d, dataset_name))

  full_template = copy.copy(multicrab)
  if "EMEnriched" in dataset_path:
    full_template.template += multicrab_semie
  elif "MuEnriched" in dataset_path:
    full_template.template += multicrab_semimu
  else:
    full_template.template += multicrab_semie
    full_template.template += multicrab_semimu

  print("Creating config file for '%s'" % (dataset_path))
  print("\tName: %s" % dataset_name)
  print("\tOutput directory (semi-mu): %s" % output_dir_semimu)
  print("\tOutput directory (semi-e): %s" % output_dir_semie)
  print("")

  f = open(output_file, "w")
  f.write(full_template.substitute(ui_working_dir=ui_working_dir, dataset=dataset_path, remote_dir_semie=output_dir_semie, remote_dir_semimu=output_dir_semimu, name=dataset_name, email=email))
  f.close()

  if options.run:
    cmd = "multicrab -create -submit -cfg %s" % (output_file)
    os.system(cmd)
