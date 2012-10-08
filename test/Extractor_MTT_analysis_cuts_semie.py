import FWCore.ParameterSet.Config as cms

analysisSettings = cms.untracked.vstring(
    "VTX_Ndof_Min             4",
    
    "MET_Pt_Min               20",
    
    "MU_Pt_min_loose          10",
    "MU_Eta_max_loose         2.5",
    "MU_Iso_min               0.125",
    "MU_Pt_min                20",
    "MU_Eta_max               2.4",
    "MU_normChi2_max          10",
    "MU_nValTrackHits_min     10",
    "MU_nMatches_min          1",
    "MU_nValPixHits_min       1",
    "MU_dB_min                0.02",
    "MU_ePt_min               15",
    "MU_eEta_max              2.5",
    "MU_eEtaW_min             1.4442",
    "MU_eEtaW_max             1.5560",
    "MU_eIso_min              0.2",

    "ELE_Iso_min              0.1",
    "ELE_Pt_min               30",
    "ELE_Eta_max              2.5",
    "ELE_Zmass                91",
    "ELE_Zwin                 15",
    "ELE_dB_min               0.02",
    
    "JET_Pt_min               30",
    "JET_Eta_max              2.4",

    "JET_btag_CSVL_min        0.244",
    "JET_btag_CSVM_min        0.679",
    "JET_btag_CSVT_min        0.898",
    "JET_btag_TCHPT_min       3.41",

    "W_mass                   80.399",
    "Top_mass                 172.0",
    "W_mass_err               10",
    "Top_mass_err             15.2",
    "b_mass                   4.67",

    "doSemiMu                 0",
    "doSyst                   0",
    "systvalue                1",
    "doUseBTaginChi2          1",
    "doChoiceWKF              0",

    #Chi2
    "chi2_hadronic_top_mass           176.583",
    "chi2_leptonic_top_mass_semimu    170.884",
    "chi2_leptonic_top_mass_semie     171.265",
    "chi2_hadronic_w_mass             84.0483",
    "chi2_pt_ttbar_system             0",
    "chi2_ht_frac                     1",

    "chi2_sigma_hadronic_top_mass           20.398",
    "chi2_sigma_leptonic_top_mass_semimu    20.4832",
    "chi2_sigma_leptonic_top_mass_semie     20.7838",
    "chi2_sigma_hadronic_w_mass             10.136",
    "chi2_sigma_pt_ttbar_system             56.93",
    "chi2_sigma_ht_frac                     0.151",

    "trigger                  ^HLT_Ele25_CaloIdVT_CaloIso(VL|T)_TrkId(VL|T)_TrkIsoT_TriCentralPF(NoPU)?Jet30_v[0-9]{0,2}$"
    )
