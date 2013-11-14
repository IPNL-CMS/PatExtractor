#! /usr/bin/env python

import os, copy, datetime, pwd, re

from optparse import OptionParser
parser = OptionParser()
parser.add_option("", "--create-cfg", action="store_true", dest="create_cfg", default=False, help="create config files for crab")
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
parser.add_option("", "--status", action="store_true", dest="status", default=False, help="run crab -status")
parser.add_option("", "--get", action="store_true", dest="get", default=False, help="run crab -get")
parser.add_option("", "--resubmit", action="store_true", dest="resubmit", default=False, help="run crab -resubmit bad")
parser.add_option("", "--submit", action="store_true", dest="submit", default=False, help="run crab -submit all")
(options, args) = parser.parse_args()

if options.run:
  options.create_cfg = True

datasets = [
    # Single Top
    ["/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/jruizalv-T_tW-channel_START53_V7A_12Nov13-v1-37c7db7f214621ff15b94bc076828bf1/USER", "T-tw"],
    ["/T_s-channel_TuneZ2star_8TeV-powheg-tauola/jruizalv-T_s-channel_START53_V7A_12Nov13-v1-37c7db7f214621ff15b94bc076828bf1/USER", "T-s"] 
    
    ]

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
CMSSW.total_number_of_events = ${events}
"""
)

multicrab_singleTprime = r"""
[${name}_singleTprime]
CMSSW.pset = Extractor_SingleTprime_MC.py
USER.user_remote_dir = ${remote_dir}
"""

for dataset in datasets:

  dataset_name = dataset[1]
  dataset_path = dataset[0]
  dataset_size = -1
  if len(dataset) > 2:
    dataset_size = dataset[2]
  #dataset_globaltag = re.search('START\d{0,2}_V\d[A-Z]?', dataset_path).group(0)

  #publish_name = "%s_%s_%s-v%d" % (dataset_name, dataset_globaltag, d, version)
  output_file = "%s_%s_multicrab_MC.cfg" % (dataset_name, d)
  ui_working_dir = ("%s_multicrab_MC") % (dataset_name)

  if options.create_cfg:
    output_dir_singleTprime = ("jruizalv/Extracted_MC/%s/%s" % (d, dataset_name))

    full_template = copy.copy(multicrab)
    full_template.template += multicrab_singleTprime

    print("Creating config file for '%s'" % (dataset_path))
    print("\tName: %s" % dataset_name)
    print("\tOutput directory (singleTprime): %s" % output_dir_singleTprime)
    print("")

    f = open(output_file, "w")
    f.write(full_template.substitute(ui_working_dir=ui_working_dir, dataset=dataset_path, remote_dir=output_dir_singleTprime, name=dataset_name, email=email, events=dataset_size))
    f.close()

  if options.run:
    cmd = "multicrab -create -submit -cfg %s" % (output_file)
    os.system(cmd)

  if options.status:
    cmd = "multicrab -status -c %s" % (ui_working_dir)
    os.system(cmd)

  if options.get:
    cmd = "multicrab -get -c %s" % (ui_working_dir)
    os.system(cmd)

  if options.resubmit:
    cmd = "multicrab -resubmit bad -c %s" % (ui_working_dir)
    os.system(cmd)

  if options.submit:
    cmd = "multicrab -submit all -c %s" % (ui_working_dir)
    os.system(cmd) 
