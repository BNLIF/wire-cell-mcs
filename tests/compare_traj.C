#include <vector>
void compare_traj(){
  TChain *T_true = new TChain("T","T");
  T_true->Add("mcs-tracks.root");
  std::vector<double> *x = new std::vector<double>;
  std::vector<double> *y = new std::vector<double>;
  std::vector<double> *z = new std::vector<double>;
  std::vector<double> *Q = new std::vector<double>;  
  Int_t N;
  T_true->SetBranchAddress("N",&N);
  T_true->SetBranchAddress("x",&x);
  T_true->SetBranchAddress("y",&y);
  T_true->SetBranchAddress("z",&z);
  T_true->SetBranchAddress("Q",&Q);
  T_true->GetEntry(0);
  // 
  //std::cout << N << " " << x->size() << std::endl;

  TGraph *g1_xy = new TGraph();
  TGraph *g1_xz = new TGraph();
  TGraph *g1_yz = new TGraph();

  for (size_t i=0;i!=x->size();i++){
    g1_xy->SetPoint(i,x->at(i),y->at(i));
    g1_xz->SetPoint(i,x->at(i),z->at(i));
    g1_yz->SetPoint(i,y->at(i),z->at(i));
  }
  g1_xy->SetLineColor(1);
  g1_xz->SetLineColor(1);
  g1_yz->SetLineColor(1);

  g1_xy->SetLineWidth(2);
  g1_xz->SetLineWidth(2);
  g1_yz->SetLineWidth(2);


  TChain *T_rec = new TChain("T_rec_charge","T_rec_charge");
  T_rec->Add("tracking_0_0_0.root");
  Double_t x1,y1,z1;
  T_rec->SetBranchAddress("x",&x1);
  T_rec->SetBranchAddress("y",&y1);
  T_rec->SetBranchAddress("z",&z1);
  TGraph *g2_xy = new TGraph();
  TGraph *g2_xz = new TGraph();
  TGraph *g2_yz = new TGraph();
  for (int i=0;i!=T_rec->GetEntries();i++){
    T_rec->GetEntry(i);

    // binning effect 1 us later from the binned slice effect, rebinned 4 ...
    // speed of imaging 1.101 mm / us
    // speed of simulation 1.098 mm/us
    // there is a potential 1 us offset at SP
    // -0.6 cm is the distance difference between Y and U planes
    x1 = (x1+0.1101*2)/1.1009999*1.098-0.6;//+ 4*0.1101; 
    
    g2_xy->SetPoint(i,x1,y1);
    g2_xz->SetPoint(i,x1,z1);
    g2_yz->SetPoint(i,y1,z1);
  }
  g2_xy->SetLineColor(2);
  g2_xz->SetLineColor(2);
  g2_yz->SetLineColor(2);

  g2_xy->SetLineWidth(1);
  g2_xz->SetLineWidth(1);
  g2_yz->SetLineWidth(1);

 
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->Divide(3,1);
  c1->cd(1);
  g1_xy->Draw("AL");
  g1_xy->GetXaxis()->SetTitle("X (cm)");
  g1_xy->GetYaxis()->SetTitle("Y (cm)");
  g2_xy->Draw("*Lsame");
  //  T_true->Draw("y:z");
  c1->cd(2);
  g1_xz->Draw("AL");
  g1_xz->GetXaxis()->SetTitle("X (cm)");
  g1_xz->GetYaxis()->SetTitle("Z (cm)");
  g2_xz->Draw("*Lsame");
  //T_rec->Draw("y:z");
  c1->cd(3);
  g1_yz->Draw("AL");
  g1_yz->GetXaxis()->SetTitle("Y (cm)");
  g1_yz->GetYaxis()->SetTitle("Z (cm)");
  g2_yz->Draw("*Lsame");

  g2_xy->SetMarkerStyle(23);
  g2_xy->SetMarkerColor(2);
  g2_xy->SetMarkerSize(0.6);

  g2_xz->SetMarkerStyle(23);
  g2_xz->SetMarkerColor(2);
  g2_xz->SetMarkerSize(0.6);

  g2_yz->SetMarkerStyle(23);
  g2_yz->SetMarkerColor(2);
  g2_yz->SetMarkerSize(0.6);
  
  TLegend *le1 = new TLegend(0.5,0.15,0.89,0.5);
  le1->SetFillColor(10);
  le1->AddEntry(g1_yz,"True","l");
  le1->AddEntry(g2_yz,"Reco","l");
  le1->Draw();
}
