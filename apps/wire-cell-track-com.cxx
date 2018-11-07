#include <iostream>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"

#include "TGraph.h"
#include "TGraph2D.h"



//using namespace WireCell;

int main(int argc, char* argv[])
{
  if (argc < 2){
    std::cerr << "usage: wire-cell-track-com -a[truth.root] -b[reco.root] -t[reco_treename] -n[truth_treename] -o[out.root]" << std::endl;
    return 1;
  }

  TString truth_filename = "mcs-tracks.root"; 
  TString reco_filename = "tracking_0_0_0.root"; 
  TString reco_treename = "T_rec_charge";
  TString truth_treename = "T";
  TString out_filename = "track_com.root";

  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 'a':
       truth_filename = &argv[i][2];
       break;
     case 'b':
       reco_filename = &argv[i][2];
       break;
     case 't':
       reco_treename = &argv[i][2];
       break;
     case 'n':
       truth_treename = &argv[i][2]; 
       break;
     case 'o':
       out_filename = &argv[i][2];
       break;
     }
  }

  std::cout << truth_filename << " " << reco_filename << " " << truth_treename << " " << reco_treename << std::endl;

  TFile *file = new TFile(out_filename,"RECREATE");
  

  
  TChain *T_true = new TChain(truth_treename,truth_treename);
  T_true->Add(truth_filename);
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

  //std::cout << N << std::endl;

  TGraph *g1_xy = new TGraph();
  TGraph *g1_xz = new TGraph();
  TGraph *g1_yz = new TGraph();

  TGraph2D *g1 = new TGraph2D();

  for (size_t i=0;i!=x->size();i++){
    g1_xy->SetPoint(i,x->at(i),y->at(i));
    g1_xz->SetPoint(i,x->at(i),z->at(i));
    g1_yz->SetPoint(i,y->at(i),z->at(i));

    g1->SetPoint(i,x->at(i),y->at(i),z->at(i));
  }
  g1_xy->SetLineColor(1);
  g1_xz->SetLineColor(1);
  g1_yz->SetLineColor(1);

  g1_xy->SetLineWidth(2);
  g1_xz->SetLineWidth(2);
  g1_yz->SetLineWidth(2);

  g1->SetLineColor(1);
  g1->SetLineWidth(2);

  TChain *T_rec = new TChain(reco_treename,reco_treename);
  T_rec->Add(reco_filename);
  Double_t x1,y1,z1;
  T_rec->SetBranchAddress("x",&x1);
  T_rec->SetBranchAddress("y",&y1);
  T_rec->SetBranchAddress("z",&z1);
  TGraph *g2_xy = new TGraph();
  TGraph *g2_xz = new TGraph();
  TGraph *g2_yz = new TGraph();

  TGraph2D *g2 = new TGraph2D();
  
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

    g2->SetPoint(i,x1,y1,z1);
    
  }

  g2_xy->SetLineColor(2);
  g2_xz->SetLineColor(2);
  g2_yz->SetLineColor(2);

  g2_xy->SetLineWidth(1);
  g2_xz->SetLineWidth(1);
  g2_yz->SetLineWidth(1);

  g2->SetLineColor(2);
  g2->SetLineWidth(2);

  //  g1_xy->GetXaxis()->SetTitle("X (cm)");
  //   g1_xy->GetYaxis()->SetTitle("Y (cm)");
  
  //    g1_xz->GetXaxis()->SetTitle("X (cm)");
  //   g1_xz->GetYaxis()->SetTitle("Z (cm)");
  
  // g1_yz->GetXaxis()->SetTitle("Y (cm)");
  //   g1_yz->GetYaxis()->SetTitle("Z (cm)");
  
  
  g2_xy->SetMarkerStyle(23);
  g2_xy->SetMarkerColor(2);
  g2_xy->SetMarkerSize(0.6);

  g2_xz->SetMarkerStyle(23);
  g2_xz->SetMarkerColor(2);
  g2_xz->SetMarkerSize(0.6);

  g2_yz->SetMarkerStyle(23);
  g2_yz->SetMarkerColor(2);
  g2_yz->SetMarkerSize(0.6);

  
  g1_xy->Write("gt_xy");
  g1_xz->Write("gt_xz");
  g1_yz->Write("gt_yz");

  g2_xy->Write("gr_xy");
  g2_xz->Write("gr_xz");
  g2_yz->Write("gr_yz");

  
  g1->Write("gt_3D");
  g2->Write("gr_3D");

  T_true->CloneTree(-1,"fast");
  T_rec->CloneTree(-1,"fast");
  
  
  file->Write();
  file->Close();
  
  
}
