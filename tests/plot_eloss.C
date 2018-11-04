void plot_eloss(){
  TFile *file = new TFile("eloss.root");
  
  Double_t T, L, dEodx;
  Double_t L_f;
  TTree *T_p = (TTree*)file->Get("T_p");
  T_p->SetBranchAddress("L",&L);
  T_p->SetBranchAddress("T",&T);
  T_p->SetBranchAddress("dEodx",&dEodx);

  TTree *T_pi = (TTree*)file->Get("T_pi");
  T_pi->SetBranchAddress("L",&L);
  T_pi->SetBranchAddress("T",&T);
  T_pi->SetBranchAddress("dEodx",&dEodx);

  TTree *T_k = (TTree*)file->Get("T_k");
  T_k->SetBranchAddress("L",&L);
  T_k->SetBranchAddress("T",&T);
  T_k->SetBranchAddress("dEodx",&dEodx);

  TTree *T_mu = (TTree*)file->Get("T_mu");
  T_mu->SetBranchAddress("L",&L);
  T_mu->SetBranchAddress("T",&T);
  T_mu->SetBranchAddress("dEodx",&dEodx);

  TTree *T_ele = (TTree*)file->Get("T_ele");
  T_ele->SetBranchAddress("L",&L);
  T_ele->SetBranchAddress("T",&T);
  T_ele->SetBranchAddress("dEodx",&dEodx);

  // proton
  T_p->GetEntry(T_p->GetEntries()-1);
  L_f = L;
  TGraph *gp = new TGraph();
  for (int i=0;i!=T_p->GetEntries()-1;i++){
    T_p->GetEntry(i);
    gp->SetPoint(i,L_f-L,dEodx);
  }

  T_pi->GetEntry(T_pi->GetEntries()-1);
  L_f = L;
  TGraph *gpi = new TGraph();
  for (int i=0;i!=T_pi->GetEntries()-1;i++){
    T_pi->GetEntry(i);
    gpi->SetPoint(i,L_f-L,dEodx);
  }

  T_k->GetEntry(T_k->GetEntries()-1);
  L_f = L;
  TGraph *gk = new TGraph();
  for (int i=0;i!=T_k->GetEntries()-1;i++){
    T_k->GetEntry(i);
    gk->SetPoint(i,L_f-L,dEodx);
  }

  T_mu->GetEntry(T_mu->GetEntries()-1);
  L_f = L;
  TGraph *gmu = new TGraph();
  for (int i=0;i!=T_mu->GetEntries()-1;i++){
    T_mu->GetEntry(i);
    gmu->SetPoint(i,L_f-L,dEodx);
  }

  T_ele->GetEntry(T_ele->GetEntries()-1);
  L_f = L;
  TGraph *gele = new TGraph();
  for (int i=0;i!=T_ele->GetEntries()-1;i++){
    T_ele->GetEntry(i);
    gele->SetPoint(i,L_f-L,dEodx);
  }


  gmu->Draw("AL");
  gmu->GetXaxis()->SetRangeUser(0,25);
  gmu->GetYaxis()->SetRangeUser(0,20);
  gmu->SetLineColor(1);
  gmu->SetLineWidth(2);
  gmu->GetXaxis()->SetTitle("L (cm)");
  gmu->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
  
  gp->Draw("Lsame");
  gp->SetLineColor(2);
  gp->SetLineWidth(2);
  
  
  gpi->Draw("Lsame");
  gpi->SetLineColor(4);
  gpi->SetLineWidth(2);

  
  
  gk->Draw("Lsame");
  gk->SetLineColor(6);
  gk->SetLineWidth(2);

  
  gele->Draw("Lsame");
  gele->SetLineColor(8);
  gele->SetLineWidth(2);

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetFillColor(10);
  le1->SetBorderSize(0);
  le1->AddEntry(gp,"Proton","l");
  le1->AddEntry(gk,"Kaon","l");
  le1->AddEntry(gpi,"Pion","l");
  le1->AddEntry(gmu,"Muon","l");
  le1->AddEntry(gele,"Electron","l");
  
  le1->Draw();
  
}
