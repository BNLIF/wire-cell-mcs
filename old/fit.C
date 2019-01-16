{
  file = TFile::Open("mcs-tracks.root");
  TH1F* h = (TH1F*)file->Get("h");


  TF1* f1 = new TF1("f1","[0] * x/[1] * exp(- x*x/(2.*[1]))",0,0.01);

  f1->SetParameter(0,1);
  f1->SetParameter(1,0.002*0.002);
  f1->SetParName(1,"#sigma^{2}");
  h->Fit(f1);
//  f1->Draw();
}
