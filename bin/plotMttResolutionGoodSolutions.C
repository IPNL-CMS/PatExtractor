{

  tree->Draw("(mtt-MC_mtt)>>mtt");
  tree->Draw("(mtt_after_kf-MC_mtt)>>mtt_after_kf", "", "sames");
  tree->Draw("(mtt_after_kf-MC_mtt)>>mtt_after_kf_with_cut", "(kf_proba > 0.2)", "sames");
  tree->Draw("(mtt_after_kf-MC_mtt)>>mtt_after_kf_with_cut_weighted", "((kf_proba > 0.2)) * kf_proba", "sames");

  gPad->Modified();
  gPad->Update();

  TH1F *mtt = (TH1F*) gDirectory->Get("mtt");
  mtt->SetTitle("");
  mtt->SetName("m_{t#bar{t}}");
  mtt->GetXaxis()->SetTitle("m_{t#bar{t}} - m_{t#bar{t}}^{gen} (GeV) (Matched events)");

  mtt = (TH1F*) gDirectory->Get("mtt_after_kf");
  mtt->SetName("m_{t#bar{t}} (after KF)");
  mtt->SetLineColor(kRed);
  mtt->SetFillStyle(3003);
  mtt->SetFillColor(kRed);
  TPaveStats *st = (TPaveStats*) mtt->FindObject("stats");
  st->SetY2NDC(0.75);
  st->SetY1NDC(0.59);
  st->SetTextColor(kRed);

  mtt = (TH1F*) gDirectory->Get("mtt_after_kf_with_cut");
  mtt->SetName("m_{t#bar{t}} (after KF, with cut)");
  mtt->SetLineColor(kViolet);
  mtt->SetFillStyle(3003);
  mtt->SetFillColor(kViolet);
  TPaveStats *st = (TPaveStats*) mtt->FindObject("stats");
  st->SetY2NDC(0.565);
  st->SetY1NDC(0.405);
  st->SetTextColor(kViolet);

  mtt = (TH1F*) gDirectory->Get("mtt_after_kf_with_cut_weighted");
  mtt->SetName("m_{t#bar{t}} (after KF, with cut, weighted)");
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
