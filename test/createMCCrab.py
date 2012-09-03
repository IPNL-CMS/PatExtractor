#! /usr/bin/env python

import os, copy

version = 1
datasets = {
    "/TTJets_TuneZ2star_8TeV-madgraph-tauola/sbrochet-TTJets_2012_v1-265c9c69c37a8e555f9b98fa1aae946f/USER": "TTJets",
    "/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/sbrochet-DYJetsToLL_M-50_2012_PF2PAT_v1-265c9c69c37a8e555f9b98fa1aae946f/USER": "DYJetsToLL_M-50",
    "/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/sbrochet-T_tW-channel_2012_PF2PAT_v1-265c9c69c37a8e555f9b98fa1aae946f/USER": "T_tW-channel",
    "/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/sbrochet-Tbar_t-channel_2012_PF2PAT_v1-265c9c69c37a8e555f9b98fa1aae946f/USER": "Tbar_t-channel",
    "/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/sbrochet-Tbar_tW-channel_2012_PF2PAT_v1-265c9c69c37a8e555f9b98fa1aae946f/USER": "Tbar_tW-channel",
    "/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/sbrochet-Tbar_s-channel_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "Tbar_s-channel",
    "/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/sbrochet-WJetsToLNu_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "WJetsToLNu",
    "/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_20_30_EMEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "QCD_Pt_20_30_EMEnriched",
    "/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_30_80_EMEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "QCD_Pt_30_80_EMEnriched",
    "/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_80_170_EMEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "QCD_Pt_80_170_EMEnriched",
    "/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_20_MuEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER": "QCD_Pt_20_MuEnriched"
    }

from string import Template
multicrab = Template(r"""[MULTICRAB]
cfg=crab_MC.cfg.template.ipnl

[COMMON]
USER.ui_working_dir = ${ui_working_dir}
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

for dataset, ui in datasets.items():
  name = ("%s_2012_v%d") % (ui, version)
  ui_working_dir = ("crab_%s") % (name)
  output_file = "multicrab_MC_%s.cfg" % (name)
  #output_dir_semie = ("MTT/Extracted/MC/Summer12/semie/%s_v%d" % (ui, version)).replace("/", "\\/")
  #output_dir_semimu = ("MTT/Extracted/MC/Summer12/semimu/%s_v%d" % (ui, version)).replace("/", "\\/")
  output_dir_semie = ("MTT/Extracted/MC/Summer12/semie/%s_v%d" % (ui, version))
  output_dir_semimu = ("MTT/Extracted/MC/Summer12/semimu/%s_v%d" % (ui, version))

  full_template = copy.copy(multicrab)
  if "EMEnriched" in dataset:
    full_template.template += multicrab_semie
  elif "MuEnriched" in dataset:
    full_template.template += multicrab_semimu
  else:
    full_template.template += multicrab_semie
    full_template.template += multicrab_semimu

  print "Creating config file for '%s'" % (dataset)
  #os.system("sed -e \"s/@datasetname@/%s/g\" -e \"s/@uiworkingdir@/%s/g\" -e \"s/@outputdir@/%s/g\" crab_MC.cfg.template.ipnl > %s" % (dataset.replace("/", "\\/"), ui_working_dir, output_dir, output_file))

  f = open(output_file, "w")
  f.write(full_template.substitute(ui_working_dir=ui_working_dir, dataset=dataset, remote_dir_semie=output_dir_semie, remote_dir_semimu=output_dir_semimu, name=name))
  f.close()
