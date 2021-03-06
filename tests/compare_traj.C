//#include <vector>
void compare_traj(){
  // TChain *T_true = new TChain("T","T");
  // T_true->Add("mcs-tracks.root");
  // std::vector<double> *x = new std::vector<double>;
  // std::vector<double> *y = new std::vector<double>;
  // std::vector<double> *z = new std::vector<double>;
  // std::vector<double> *Q = new std::vector<double>;  
  // Int_t N;
  // T_true->SetBranchAddress("N",&N);
  // T_true->SetBranchAddress("x",&x);
  // T_true->SetBranchAddress("y",&y);
  // T_true->SetBranchAddress("z",&z);
  // T_true->SetBranchAddress("Q",&Q);
  // T_true->GetEntry(0);
  // 
  //std::cout << N << " " << x->size() << std::endl;

  // TGraph *g1_xy = new TGraph();
  // TGraph *g1_xz = new TGraph();
  // TGraph *g1_yz = new TGraph();

  // for (size_t i=0;i!=x->size();i++){
  //   g1_xy->SetPoint(i,x->at(i),y->at(i));
  //   g1_xz->SetPoint(i,x->at(i),z->at(i));
  //   g1_yz->SetPoint(i,y->at(i),z->at(i));
  // }
  // g1_xy->SetLineColor(1);
  // g1_xz->SetLineColor(1);
  // g1_yz->SetLineColor(1);

  // g1_xy->SetLineWidth(2);
  // g1_xz->SetLineWidth(2);
  // g1_yz->SetLineWidth(2);

  TFile *file = new TFile("track_com.root");
  TGraph *g1_xy = (TGraph*)file->Get("gt_xy");
  TGraph *g1_xz = (TGraph*)file->Get("gt_xz");
  TGraph *g1_yz = (TGraph*)file->Get("gt_yz");
  TGraph2D *g1 = (TGraph2D*)file->Get("gt_3D");
  
  TGraph *g2_xy = (TGraph*)file->Get("gr_xy");
  TGraph *g2_xz = (TGraph*)file->Get("gr_xz");
  TGraph *g2_yz = (TGraph*)file->Get("gr_yz");
  TGraph2D *g2 = (TGraph2D*)file->Get("gr_3D");

  TTree *t1 = (TTree*)file->Get("T");
  double max_dis=0; // maximum distance along the track 
  double beg_dis; // distance at the beginning ... 
  double end_dis; // distance at the end ...
  double total_dis2=0; // total distance^2
  double total_L = 0;
  int Npoints = 0;
  double total_dtheta = 0;
  double max_dtheta = 0;
  
  t1->SetBranchAddress("max_dis",&max_dis);
  t1->SetBranchAddress("beg_dis",&beg_dis);
  t1->SetBranchAddress("end_dis",&end_dis);
  t1->SetBranchAddress("total_dis2",&total_dis2);
  t1->SetBranchAddress("total_L",&total_L);
  t1->SetBranchAddress("N",&Npoints);
  t1->Branch("total_dtheta",&total_dtheta);
  t1->Branch("max_dtheta",&max_dtheta);
  
  std::vector<double> *x2 = new std::vector<double>;
  std::vector<double> *y2 = new std::vector<double>;
  std::vector<double> *z2 = new std::vector<double>;
  std::vector<double> *dQ_rec = new std::vector<double>;
  std::vector<double> *dQ_tru = new std::vector<double>;
  std::vector<double> *dx = new std::vector<double>;
  
  std::vector<double> *x2_pair = new std::vector<double>;
  std::vector<double> *y2_pair = new std::vector<double>;
  std::vector<double> *z2_pair = new std::vector<double>;

  std::vector<double> *dis = new std::vector<double>;
  std::vector<double> *L = new std::vector<double>;
  std::vector<double> *dtheta = new std::vector<double>;
  
  t1->SetBranchAddress("x",&x2);
  t1->SetBranchAddress("y",&y2);
  t1->SetBranchAddress("z",&z2);
  t1->SetBranchAddress("dQ_rec",&dQ_rec);
  t1->SetBranchAddress("dQ_tru",&dQ_tru);
  t1->SetBranchAddress("dx",&dx);
  
  t1->SetBranchAddress("x_pair",&x2_pair);
  t1->SetBranchAddress("y_pair",&y2_pair);
  t1->SetBranchAddress("z_pair",&z2_pair);

  t1->SetBranchAddress("dis",&dis);
  t1->SetBranchAddress("L",&L);
  t1->SetBranchAddress("dtheta",&dtheta);


  t1->GetEntry(0);
  // //TChain *T_rec = new TChain("T_rec_charge","T_rec_charge");
  // TChain *T_rec = new TChain("T_rec","T_rec");
  // T_rec->Add("tracking_0_0_0.root");
  // Double_t x1,y1,z1;
  // T_rec->SetBranchAddress("x",&x1);
  // T_rec->SetBranchAddress("y",&y1);
  // T_rec->SetBranchAddress("z",&z1);
  // TGraph *g2_xy = new TGraph();
  // TGraph *g2_xz = new TGraph();
  // TGraph *g2_yz = new TGraph();
  // for (int i=0;i!=T_rec->GetEntries();i++){
  //   T_rec->GetEntry(i);

  //   // binning effect 1 us later from the binned slice effect, rebinned 4 ...
  //   // speed of imaging 1.101 mm / us
  //   // speed of simulation 1.098 mm/us
  //   // there is a potential 1 us offset at SP
  //   // -0.6 cm is the distance difference between Y and U planes
  //   x1 = (x1+0.1101*2)/1.1009999*1.098-0.6;//+ 4*0.1101; 
    
  //   g2_xy->SetPoint(i,x1,y1);
  //   g2_xz->SetPoint(i,x1,z1);
  //   g2_yz->SetPoint(i,y1,z1);
  // }
  // g2_xy->SetLineColor(2);
  // g2_xz->SetLineColor(2);
  // g2_yz->SetLineColor(2);

  // g2_xy->SetLineWidth(1);
  // g2_xz->SetLineWidth(1);
  // g2_yz->SetLineWidth(1);

 
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  c1->Divide(4,2);
  c1->cd(1);
  g1->Draw("AP");
  g2->Draw("Psame");

  g1->SetMarkerColor(1);
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.5);

  g2->SetMarkerColor(2);
  g2->SetMarkerStyle(21);
  g2->SetMarkerSize(0.5);

  g1->GetXaxis()->SetTitle("X (cm)");
  g1->GetYaxis()->SetTitle("Y (cm)");
  g1->GetZaxis()->SetTitle("Z (cm)");
  g1->SetTitle("3D trajectory");

  TPolyMarker3D *gm = new TPolyMarker3D();
  gm->SetPoint(0,x2->front(),y2->front(),z2->front());
  gm->SetMarkerColor(4);
  gm->SetMarkerStyle(29);
  gm->SetMarkerSize(2);
  gm->Draw("*same");
  
  c1->cd(2);
  g1_xy->Draw("AL");
  g1_xy->GetXaxis()->SetTitle("X (cm)");
  g1_xy->GetYaxis()->SetTitle("Y (cm)");
  g2_xy->Draw("*Lsame");

  TMarker *p1 = new TMarker(x2->front(), y2->front(), 29);
  p1->Draw("same");
  p1->SetMarkerColor(4);
  p1->SetMarkerSize(2);
  //  T_true->Draw("y:z");
  c1->cd(5);
  g1_xz->Draw("AL");
  g1_xz->GetXaxis()->SetTitle("X (cm)");
  g1_xz->GetYaxis()->SetTitle("Z (cm)");
  g2_xz->Draw("*Lsame");

  TMarker *p2 = new TMarker(x2->front(), z2->front(), 29);
  p2->Draw("same");
  p2->SetMarkerColor(4);
  p2->SetMarkerSize(2);
  //T_rec->Draw("y:z");
  c1->cd(6);
  g1_yz->Draw("AL");
  g1_yz->GetXaxis()->SetTitle("Y (cm)");
  g1_yz->GetYaxis()->SetTitle("Z (cm)");
  g2_yz->Draw("*Lsame");

  TMarker *p3 = new TMarker(y2->front(), z2->front(), 29);
  p3->Draw("same");
  p3->SetMarkerColor(4);
  p3->SetMarkerSize(2);

  
  
  TLegend *le1 = new TLegend(0.5,0.15,0.89,0.5);
  le1->SetFillColor(10);
  le1->AddEntry(g1_yz,"True","l");
  le1->AddEntry(g2_yz,"Reco","l");
  le1->Draw();


  c1->cd(3);
  TGraph *gd = new TGraph();
  for (size_t i=0;i!=dis->size();i++){
    gd->SetPoint(i,L->at(i),dis->at(i));
  }
  gd->Draw("APL");
  gd->SetMarkerStyle(20);
  gd->GetXaxis()->SetTitle("Track Length (cm)");
  gd->GetYaxis()->SetTitle("distance (cm)");
  
  c1->cd(4);
  //T->Draw("dtheta:L");
  TGraph *gt = new TGraph();
  for (size_t i=0;i!=dtheta->size();i++){
    gt->SetPoint(i,L->at(i),dtheta->at(i)/3.1415926*180.);
  }
  gt->Draw("APL");
  gt->SetMarkerStyle(20);
  gt->GetXaxis()->SetTitle("Track Length (cm)");
  gt->GetYaxis()->SetTitle("angle (deg)");

  c1->cd(7);
  TGraph *g_dQdx = new TGraph();
  TGraph *g_dQdx1 = new TGraph();
  for (size_t i=0;i!=dx->size();i++){
    g_dQdx->SetPoint(i,L->at(i),dQ_rec->at(i)/dx->at(i));
    g_dQdx1->SetPoint(i,L->at(i),dQ_tru->at(i)/dx->at(i));
  }
  g_dQdx->Draw("APL");
  g_dQdx1->Draw("Psame");
  g_dQdx1->SetMarkerStyle(25);
  g_dQdx1->SetMarkerColor(1);
  g_dQdx->SetMarkerColor(2);
  g_dQdx->SetLineColor(2);
  g_dQdx->SetMarkerStyle(20);
  g_dQdx->GetXaxis()->SetTitle("Track Length (cm)");
  g_dQdx->GetYaxis()->SetTitle("dQ/dx (e/cm)");
  
}
