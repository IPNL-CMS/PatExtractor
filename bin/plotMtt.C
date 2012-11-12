{

  tree->Draw("mtt>>mtt");
  tree->Draw("mtt_after_kf>>mtt_after_kf", "kf_converged", "sames");
  tree->Draw("mtt_after_kf>>mtt_after_kf_with_cut", "kf_converged && (kf_proba > 0.2)", "sames");
  tree->Draw("mtt_after_kf>>mtt_after_kf_with_cut_weighted", "(kf_converged && (kf_proba > 0.2)) * kf_proba", "sames");

  gPad->Modified();
  gPad->Update();

  TH1F *mtt = (TH1F*) gDirectory->Get("mtt");
  mtt->GetXaxis()->SetTitle("Mass (GeV)");

  mtt = (TH1F*) gDirectory->Get("mtt_after_kf");
  mtt->SetLineColor(kRed);
  mtt->SetFillStyle(3003);
  mtt->SetFillColor(kRed);
  TPaveStats *st = (TPaveStats*) mtt->FindObject("stats");
  st->SetY2NDC(0.75);
  st->SetY1NDC(0.59);
  st->SetTextColor(kRed);

  mtt = (TH1F*) gDirectory->Get("mtt_after_kf_with_cut");
  mtt->SetLineColor(kViolet);
  mtt->SetFillStyle(3003);
  mtt->SetFillColor(kViolet);
  TPaveStats *st = (TPaveStats*) mtt->FindObject("stats");
  st->SetY2NDC(0.565);
  st->SetY1NDC(0.405);
  st->SetTextColor(kViolet);

  mtt = (TH1F*) gDirectory->Get("mtt_after_kf_with_cut_weighted");
  mtt->SetLineColor(kGreen + 2);
  mtt->SetFillStyle(3003);
  mtt->SetFillColor(kGreen + 2);
  TPaveStats *st = (TPaveStats*) mtt->FindObject("stats");
  st->SetY2NDC(0.38);
  st->SetY1NDC(0.22);
  st->SetTextColor(kGreen + 2);

  gPad->Modified();
  gPad->Update();
}
