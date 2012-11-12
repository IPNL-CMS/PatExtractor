void plot(TTree* tree, TString a, TString b, float cut, TString outputName) {

  TString cutA = TString::Format("mtt_after_kf > %.02f", cut);
  TString cutB = TString::Format("mtt_after_kf < %.02f", cut);

  TString drawA = TString::Format("%s - %s>>inside_peak", a.Data(), b.Data());
  TString drawB = TString::Format("%s - %s>>outside_peak", a.Data(), b.Data());

  tree->Draw(drawA, cutA);
  tree->Draw(drawB, cutB, "sames");

  gPad->Modified();
  gPad->Update();

  TH1F *h1 = (TH1F*) gDirectory->Get("inside_peak");

  TH1F *h2 = (TH1F*) gDirectory->Get("outside_peak");
  h2->SetLineColor(kRed);
  TPaveStats *st = (TPaveStats*) h2->FindObject("stats");
  st->SetY2NDC(0.75);
  st->SetY1NDC(0.59);
  st->SetTextColor(kRed);

  gPad->Modified();
  gPad->Update();

  c1->Print(outputName);

  delete h1;
  delete h2;
}

void plotDeltaGoodSolutions(float cut = 600) {

  plot(tree, "firstJet_pt", "firstJet_after_kf_pt", cut, "kinfit_firstJet_pt_delta.pdf");
  plot(tree, "firstJet_eta", "firstJet_after_kf_eta", cut, "kinfit_firstJet_eta_delta.pdf");
  plot(tree, "firstJet_phi", "firstJet_after_kf_phi", cut, "kinfit_firstJet_phi_delta.pdf");

  plot(tree, "secondJet_pt", "secondJet_after_kf_pt", cut, "kinfit_secondJet_pt_delta.pdf");
  plot(tree, "secondJet_eta", "secondJet_after_kf_eta", cut, "kinfit_secondJet_eta_delta.pdf");
  plot(tree, "secondJet_phi", "secondJet_after_kf_phi", cut, "kinfit_secondJet_phi_delta.pdf");

  plot(tree, "hadronicB_pt", "hadronicB_after_kf_pt", cut, "kinfit_hadronicB_pt_delta.pdf");
  plot(tree, "hadronicB_eta", "hadronicB_after_kf_eta", cut, "kinfit_hadronicB_eta_delta.pdf");
  plot(tree, "hadronicB_phi", "hadronicB_after_kf_phi", cut, "kinfit_hadronicB_phi_delta.pdf");

  plot(tree, "leptonicB_pt", "leptonicB_after_kf_pt", cut, "kinfit_leptonicB_pt_delta.pdf");
  plot(tree, "leptonicB_eta", "leptonicB_after_kf_eta", cut, "kinfit_leptonicB_eta_delta.pdf");
  plot(tree, "leptonicB_phi", "leptonicB_after_kf_phi", cut, "kinfit_leptonicB_phi_delta.pdf");

  /*
  plot(tree, "lepton_pt", "lepton_after_kf_pt", cut, "kinfit_lepton_pt_delta.pdf");
  plot(tree, "lepton_eta", "lepton_after_kf_eta", cut, "kinfit_lepton_eta_delta.pdf");
  plot(tree, "lepton_phi", "lepton_after_kf_phi", cut, "kinfit_lepton_phi_delta.pdf");
  */

  plot(tree, "neutrino_pt", "neutrino_after_kf_pt", cut, "kinfit_neutrino_pt_delta.pdf");
  plot(tree, "neutrino_eta", "neutrino_after_kf_eta", cut, "kinfit_neutrino_eta_delta.pdf");
  plot(tree, "neutrino_phi", "neutrino_after_kf_phi", cut, "kinfit_neutrino_phi_delta.pdf");
  plot(tree, "neutrino_pz", "neutrino_after_kf_pz", cut, "kinfit_neutrino_pz_delta.pdf");
}
