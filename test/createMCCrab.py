#! /usr/bin/env python

import os, copy, datetime

datasets = {
    "/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/jruizalv-T_tW-channel_START53_V7A_12Nov13-v1-37c7db7f214621ff15b94bc076828bf1/USER" : "T-tW",
    "/T_s-channel_TuneZ2star_8TeV-powheg-tauola/jruizalv-T_s-channel_START53_V7A_12Nov13-v1-37c7db7f214621ff15b94bc076828bf1/USER" : "T-s"
    }

d = datetime.datetime.now().strftime("%d%b")

from string import Template
multicrab = Template(r"""[MULTICRAB]
cfg=crab_MC.cfg.template.ipnl

[COMMON]
USER.ui_working_dir = ${ui_working_dir}
CMSSW.datasetpath = ${dataset}
"""
)

multicrab_singleTprime = r"""
[${name}]
CMSSW.pset = Extractor_SingleTprime_MC.py
USER.user_remote_dir = ${remote_dir}
"""

for dataset, ui in datasets.items():
  name = ("%s_2013_%s") % (ui, d)
  ui_working_dir = ("%s_crab") % (name)
  output_file = "%s_multicrab_MC.cfg" % (name)

  output_dir_singleTprime = ("jruizalv/Extracted_MC/%s/%s" % (d, ui))

  full_template = copy.copy(multicrab)
  full_template.template += multicrab_singleTprime
  
  print "Creating config file for '%s'" % (dataset)
  #os.system("sed -e \"s/@datasetname@/%s/g\" -e \"s/@uiworkingdir@/%s/g\" -e \"s/@outputdir@/%s/g\" crab_MC.cfg.template.ipnl > %s" % (dataset.replace("/", "\\/"), ui_working_dir, output_dir, output_file))

  f = open(output_file, "w")
  f.write(full_template.substitute(ui_working_dir=ui_working_dir, dataset=dataset, remote_dir=output_dir_singleTprime, name=name))
  f.close()
