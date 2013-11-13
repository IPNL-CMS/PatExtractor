#########################################
#
# Base macro for launching the PatExtractor
#
# The macro is for tests
#
#########################################


import FWCore.ParameterSet.Config as cms

def readFile(file):
  return cms.untracked.string(open(file).read())

def createExtractorProcess(isMC, isSemiMu, useShiftCorrectedMET, globalTag):
  process = cms.Process("PATextractor2")


  #########################################
  #
  # Main configuration statements
  #
  #########################################

  process.load('Configuration/StandardSequences/Services_cff')
  process.load('Configuration/StandardSequences/GeometryIdeal_cff')
  process.load('Configuration/StandardSequences/MagneticField_38T_cff')
  process.load('Configuration/StandardSequences/EndOfProcess_cff')
  process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
  process.load("FWCore.MessageLogger.MessageLogger_cfi")
  process.load("Extractors.PatExtractor.PAT_extractor_cff")

  process.maxEvents = cms.untracked.PSet(
      input = cms.untracked.int32(10) #
      )

  #Global tag and data type choice
  process.GlobalTag.globaltag = '%s::All' % globalTag
  process.PATextraction.isMC  = isMC
  process.PATextraction.doMC  = isMC

  #Input PAT file to extract
  process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(
        ),
      duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' ),
      )

  #Output extracted file name
  if isMC:
    process.PATextraction.extractedRootFile = cms.string('extracted_mc.root')
  else:
    process.PATextraction.extractedRootFile = cms.string('extracted.root')

  #########################################
  #
  # PAT extractor main options statements
  #
  #########################################

  #
  # Adapt it to your needs
  #
  # If you are lost, see the example here (PART 3.2):
  # http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.PHYTuto
  #
  # Here we just extract, and don't perform any analysis

  process.PATextraction.doMuon     = True
  process.PATextraction.doElectron = True
  process.PATextraction.doJet      = True

  process.PATextraction.doMET      = True
  if useShiftCorrectedMET:
    process.PATextraction.MET_PF.input = cms.InputTag("patMETsShiftCorrPFlow")
  else:
    process.PATextraction.MET_PF.input  = cms.InputTag("patMETsPFlow")

  process.PATextraction.doVertex   = True
  process.PATextraction.vtx_tag    = cms.InputTag( "goodOfflinePrimaryVertices" )
  process.PATextraction.doHLT      = True

  if not isMC:
    if isSemiMu:
      process.PATextraction.triggersXML = readFile("triggers_mu.xml")
    else:
      process.PATextraction.triggersXML = readFile("triggers_e.xml")

  process.PATextraction.doMtt      = True

  # Jets correction : needs a valid global tags, or an external DB where JEC are stored
  process.PATextraction.jet_PF.redoJetCorrection = False

  if isMC:
    process.PATextraction.jet_PF.jetCorrectorLabel = "ak5PFchsL1FastL2L3"
  else:
    process.PATextraction.jet_PF.jetCorrectorLabel = "ak5PFchsL1FastL2L3Residual"

  process.PATextraction.jet_PF.doJER = False # Disable automatically on data

  # JER systematics:
  # Use -1 for 1-sigma down, 0 for nominal correction, and 1 for 1-sigma up
  process.PATextraction.jet_PF.jerSign = 0

  # JES systematics:
  # Use -1 for 1-sigma down, 0 for nominal correction, and 1 for 1-sigma up
  process.PATextraction.jet_PF.jesSign = 0
  process.PATextraction.jet_PF.jes_uncertainties_file = cms.untracked.string("Extractors/PatExtractor/data/START53_V23_Uncertainty_AK5PFchs.txt")

  process.PATextraction.MET_PF.redoMetPhiCorrection   = False
  process.PATextraction.MET_PF.redoMetTypeICorrection = False # Automatically true if redoJetCorrection is True

  from Extractor_MTT_ScaleFactors import loadMuonScaleFactor, loadBTagScaleFactors, loadElectronScaleFactor, loadLightJetsScaleFactor
  loadBTagScaleFactors(process)

  # Scale factors
  process.PATextraction.muon_scale_factors_tighteff_tightiso = loadMuonScaleFactor("MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl", "MuonEfficiencies_Run2012ReReco_53X.pkl", "Tight", "combRelIsoPF04dBeta<012_Tight")
  process.PATextraction.muon_scale_factors_tighteff_looseiso = loadMuonScaleFactor("MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl", "MuonEfficiencies_Run2012ReReco_53X.pkl", "Tight", "combRelIsoPF04dBeta<02_Tight")
  process.PATextraction.muon_scale_factors_looseeff_looseiso = loadMuonScaleFactor("MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl", "MuonEfficiencies_Run2012ReReco_53X.pkl", "Loose", "combRelIsoPF04dBeta<02_Loose")
  process.PATextraction.muon_scale_factors = cms.vstring("muon_scale_factors_looseeff_looseiso", "muon_scale_factors_tighteff_looseiso", "muon_scale_factors_tighteff_tightiso")

  process.PATextraction.electron_scale_factors_tighteff_tightiso = loadElectronScaleFactor("Electron_scale_factors.json", "tight")
  process.PATextraction.electron_scale_factors_looseeff_tightiso = loadElectronScaleFactor("Electron_scale_factors.json", "loose")
  process.PATextraction.electron_scale_factors = cms.vstring("electron_scale_factors_tighteff_tightiso", "electron_scale_factors_looseeff_tightiso")

  process.PATextraction.b_tagging_scale_factors_b_jets = cms.PSet(
          jet_type = cms.string("b"),
          from_globaltag = cms.bool(True),
          payload = cms.string("MUJETSWPBTAGTTBARCSVM")
          )
  process.PATextraction.b_tagging_scale_factors_c_jets = cms.PSet(
          jet_type = cms.string("c"),
          from_globaltag = cms.bool(True),
          payload = cms.string("MUJETSWPBTAGTTBARCSVM")
          )
  process.PATextraction.b_tagging_scale_factors_light_jets = cms.PSet(
          jet_type = cms.string("light"),
          from_globaltag = cms.bool(False),
          scale_factors = loadLightJetsScaleFactor()
          )
  process.PATextraction.b_tagging_scale_factors = cms.vstring("b_tagging_scale_factors_b_jets", "b_tagging_scale_factors_c_jets", "b_tagging_scale_factors_light_jets")

  # MTT analysis configuration
  process.PATextraction.plugins = cms.PSet(
      mtt_analysis = cms.PSet(
        do_semimu = cms.bool(isSemiMu),
        do_pdf_systematics = cms.untracked.bool(False),
        met = cms.PSet(
          pt_min = cms.double(20)
          ),

        muons_tight = cms.PSet(
          pt_min = cms.double(27),
          eta_max = cms.double(2.1),
          isolation_max = cms.double(0.12)
          ),

        muons_loose = cms.PSet(
          pt_min = cms.double(10),
          eta_max = cms.double(2.5),
          isolation_max = cms.double(0.20)
          ),

        electrons_tight = cms.PSet(
          pt_min = cms.double(30),
          eta_max = cms.double(2.5),
          isolation_max = cms.double(0.10)
          ),

        electrons_loose = cms.PSet(
          pt_min = cms.double(20),
          eta_max = cms.double(2.5),
          isolation_max = cms.double(0.15)
          ),

        jets = cms.PSet(
          pt_min = cms.double(30),
          eta_max = cms.double(2.4),
          btag_CSVL = cms.double(0.244),
          btag_CSVM = cms.double(0.679),
          btag_CSVT = cms.double(0.898),
          btag_TCHPT = cms.double(3.41),
          b_tagging_efficiency = cms.double(0.6915)
          ),

        chi2_sorting = cms.PSet(
          w_mass = cms.double(80.399),
          w_mass_error = cms.double(10),
          top_mass = cms.double(172),
          top_mass_error = cms.double(10),
          b_mass = cms.double(4.67),

          use_btagging = cms.bool(False),

          hadronic_top_mass = cms.double(175.164),
          leptonic_top_mass_semimu = cms.double(170.937),
          leptonic_top_mass_semie = cms.double(170.875),
          hadronic_w_mass = cms.double(84.0569),
          pt_ttbar_system = cms.double(0),
          ht_frac = cms.double(1),

          sigma_hadronic_top_mass = cms.double(17.3484),
          sigma_leptonic_top_mass_semimu = cms.double(17.3594),
          sigma_leptonic_top_mass_semie = cms.double(17.2875),
          sigma_hadronic_w_mass = cms.double(10.1166),
          sigma_pt_ttbar_system = cms.double(56.93),
          sigma_ht_frac = cms.double(0.151),

          use_pt_syst = cms.bool(False),
          use_ht_frac = cms.bool(False)
          ),

        mva = cms.PSet(
          weights = cms.string("Extractors/PatExtractor/data/ttbar_BDT.weights.xml")
          ),

        use_mva = cms.bool(True),
        use_chi2 = cms.bool(True),

        b_tagging_efficiency = cms.double(0.6915),

        # ------------------------------------------------
        # settings for the KinFitter
        # ------------------------------------------------
        maxNrIter = cms.uint32(500),
        maxDeltaS = cms.double(5e-05),
        maxF      = cms.double(0.0001),
        # ------------------------------------------------
        # select parametrisation
        # 0: EMom, 1: EtEtaPhi, 2: EtThetaPhi
        # ------------------------------------------------
        jetParametrisation = cms.uint32(1),
        lepParametrisation = cms.uint32(1),
        metParametrisation = cms.uint32(1),

        # ------------------------------------------------
        # set constraints
        # 1: Whad-mass, 2: Wlep-mass, 3: thad-mass,
        # 4: tlep-mass, 5: nu-mass, 6: equal t-masses
        # 7: sum-pt conservation
        # ------------------------------------------------
        constraints = cms.vuint32(1, 2, 3, 4),

        # ------------------------------------------------
        # set mass values used in the constraints
        # ------------------------------------------------
        mW   = cms.double(80.4),
        mTop = cms.double(173.),

        # ------------------------------------------------
        # set correction factor(s) for the jet energy resolution:
        # - (optional) eta dependence assumed to be symmetric
        #   around eta=0, i.e. parametrized in |eta|
        # - any negative value as last bin edge is read as "inf"
        # - make sure that number of entries in vector with
        #   bin edges = number of scale factors + 1
        # ------------------------------------------------
        jetEnergyResolutionScaleFactors = cms.vdouble(1.0),
        jetEnergyResolutionEtaBinning = cms.vdouble(0.0,-1.0))
      )

  #########################################
  #
  # Launch the job
  #
  #########################################


  process.p = cms.Path(process.PATextraction)
  process.MessageLogger.cerr.FwkReport.reportEvery = 1000

  return process
