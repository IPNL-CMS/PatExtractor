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
    # Z' narrow
    ["/ZPrimeToTTJets_M500GeV_W5GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_500_Narrow_START53_V7A_04Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_500_Narrow"],
    ["/ZPrimeToTTJets_M750GeV_W7p5GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_750_Narrow_START53_V7A_04Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_750_Narrow"],
    ["/ZPrimeToTTJets_M1000GeV_W10GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_1000_Narrow_START53_V7A_04Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_1000_Narrow"],
    ["/ZPrimeToTTJets_M1250GeV_W12p5GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_1250_Narrow_START53_V7A_04Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_1250_Narrow"],
    ["/ZPrimeToTTJets_M1500GeV_W15GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_1500_Narrow_START53_V7A_04Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_1500_Narrow"],
    ["/ZPrimeToTTJets_M2000GeV_W20GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_2000_Narrow_START53_V7A_04Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_2000_Narrow"],

    ["/ZPrimeToTTJets_M500GeV_W5GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_500_Narrow_ext_START53_V7C_28Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_500_Narrow_ext"],
    ["/ZPrimeToTTJets_M750GeV_W7p5GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_750_Narrow_ext_START53_V7C_03Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_750_Narrow_ext"],
    ["/ZPrimeToTTJets_M1000GeV_W10GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_1000_Narrow_ext_START53_V7C_03Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_1000_Narrow_ext"],
    ["/ZPrimeToTTJets_M1250GeV_W12p5GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_1250_Narrow_ext_START53_V7C_03Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_1250_Narrow_ext"],
    ["/ZPrimeToTTJets_M1500GeV_W15GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_1500_Narrow_ext_START53_V7C_03Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_1500_Narrow_ext", 100000],
    ["/ZPrimeToTTJets_M2000GeV_W20GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_2000_Narrow_ext_START53_V7C_03Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_2000_Narrow_ext"],

    # Z' large
    ["/ZPrimeToTTJets_M500GeV_W50GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_500_Large_START53_V7A_11Dec12-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_500_Large"],
    ["/ZPrimeToTTJets_M750GeV_W75GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_750_Large_START53_V7A_11Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_750_Large"], 
    ["/ZPrimeToTTJets_M1000GeV_W100GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_1000_Large_START53_V7A_11Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_1000_Large"], 
    ["/ZPrimeToTTJets_M1250GeV_W125GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_1250_Large_START53_V7A_11Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_1250_Large"], 
    ["/ZPrimeToTTJets_M1500GeV_W150GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_1500_Large_START53_V7A_11Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_1500_Large"],
    ["/ZPrimeToTTJets_M2000GeV_W200GeV_TuneZ2star_8TeV-madgraph-tauola/sperries-Zprime_2000_Large_START53_V7A_11Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Zprime_2000_Large"],

    ["/ZPrimeToTTJets_M750GeV_W75GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_750_Large_ext_START53_V7C_03Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_750_Large_ext"], 
    ["/ZPrimeToTTJets_M1000GeV_W100GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_1000_Large_ext_START53_V7C_03Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_1000_Large_ext"], 
    ["/ZPrimeToTTJets_M1250GeV_W125GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_1250_Large_ext_START53_V7C_03Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_1250_Large_ext"], 
    ["/ZPrimeToTTJets_M1500GeV_W150GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_1500_Large_ext_START53_V7C_03Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_1500_Large_ext"],
    ["/ZPrimeToTTJets_M2000GeV_W200GeV_TuneZ2star_8TeV_ext-madgraph-tauola/sbrochet-Zprime_2000_Large_ext_START53_V7C_03Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "Zprime_2000_Large_ext", 100000],

    # RS Gluons
    ["/RSGluonToTT_M-700_Tune4C_8TeV-pythia8/sperries-RSGluon700_START53_V7A_11Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "RSGluon700"], 
    ["/RSGluonToTT_M-1000_Tune4C_8TeV-pythia8/sperries-RSGluon1000_START53_V7A_11Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "RSGluon1000"],
    ["/RSGluonToTT_M-1200_Tune4C_8TeV-pythia8/sperries-RSGluon1200_START53_V7A_11Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "RSGluon1200"], 
    ["/RSGluonToTT_M-1500_Tune4C_8TeV-pythia8/sperries-RSGluon1500_START53_V7A_11Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "RSGluon1500"],
    ["/RSGluonToTT_M-2000_Tune4C_8TeV-pythia8/sperries-RSGluon2000_START53_V7A_11Dec12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "RSGluon2000"],


    # Single anti-top
    ["/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/chassera-Tbar_t-channel_START53_V7A_22Nov12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Tbar_t-channel"],
    ["/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/chassera-Tbar_tW-channel_START53_V7A_22Nov12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Tbar_tW-channel"],
    ["/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/chassera-Tbar_s-channel_START53_V7A_22Nov12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "Tbar_s-channel"],
    
    # Single top
    ["/T_t-channel_TuneZ2star_8TeV-powheg-tauola/chassera-T_t-channel_START53_V7A_22Nov12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "T_t-channel"],
    ["/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/chassera-T_tW-channel_START53_V7A_22Nov12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "T_tW-channel"],
    ["/T_s-channel_TuneZ2star_8TeV-powheg-tauola/chassera-T_s-channel_START53_V7A_22Nov12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "T_s-channel"],

    # TT + jets
    # ["/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/chassera-TTJets_MassiveBinDECAY_START53_V7A_22Nov12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "TTJets_MassiveBinDECAY"],
    ["/TT_CT10_TuneZ2star_8TeV-powheg-tauola/sperries-TT_powheg_START53_V7A_16Jan13-v1-bdd0c9c28c68bfd05bfd28ee5e93863c/USER", "TT_powheg"],
    
    ["/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/chassera-DYJetsToLL_M-50_START53_V7A_22Nov12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "DYJetsToLL_M-50"],
    ["/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/chassera-WJetsToLNu_START53_V7A_22Nov12-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "WJetsToLNu"],
   
    #["/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_20_30_EMEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER", "QCD_Pt_20_30_EMEnriched"],
    #["/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_30_80_EMEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER", "QCD_Pt_30_80_EMEnriched"],
    #["/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_80_170_EMEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER", "QCD_Pt_80_170_EMEnriched"],
    #["/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_20_MuEnriched_2012_PF2PAT_v1-3a57158a5a24f1281931c60ed2517d66/USER", "QCD_Pt_20_MuEnriched"],

    # Dibosons
    #["/WW_TuneZ2star_8TeV_pythia6_tauola/chassera-WWincl_START53_V7A_08Jan13-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "WW_incl"],
    #["/WZ_TuneZ2star_8TeV_pythia6_tauola/chassera-WZincl_START53_V7A_08Jan13-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "WZ_incl"],
    #["/ZZ_TuneZ2star_8TeV_pythia6_tauola/chassera-ZZincl_START53_V7A_08Jan13-v1-bd09b58f34b981e2c3ef3678b9b096ed/USER", "ZZ_incl"],
    
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

for dataset in datasets:

  dataset_name = dataset[1]
  dataset_path = dataset[0]
  dataset_size = -1
  if len(dataset) > 2:
    dataset_size = dataset[2]
  #dataset_globaltag = re.search('START\d{0,2}_V\d[A-Z]?', dataset_path).group(0)

  #publish_name = "%s_%s_%s-v%d" % (dataset_name, dataset_globaltag, d, version)
  output_file = "multicrab_MC_%s_%s.cfg" % (dataset_name, d)
  ui_working_dir = ("multicrab_MC_%s") % (dataset_name)

  if options.create_cfg:
    output_dir_semie = ("Extracted_step2/MC/Summer12/%s/semie/%s" % (d, dataset_name))
    output_dir_semimu = ("Extracted_step2/MC/Summer12/%s/semimu/%s" % (d, dataset_name))

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
    f.write(full_template.substitute(ui_working_dir=ui_working_dir, dataset=dataset_path, remote_dir_semie=output_dir_semie, remote_dir_semimu=output_dir_semimu, name=dataset_name, email=email, events=dataset_size))
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
