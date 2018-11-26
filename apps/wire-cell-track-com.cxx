
#include "WireCellData/ToyPointCloud.h"

#include <iostream>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"

#include "TGraph.h"
#include "TGraph2D.h"
#include "TVector3.h"

#include <map>

using namespace WireCell;

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

  ToyPointCloud pcloud, pcloud1;

  WireCell::PointVector ps;
  
  for (size_t i=0;i!=x->size();i++){
    g1_xy->SetPoint(i,x->at(i),y->at(i));
    g1_xz->SetPoint(i,x->at(i),z->at(i));
    g1_yz->SetPoint(i,y->at(i),z->at(i));

    Point p(x->at(i)*units::cm,y->at(i)*units::cm,z->at(i)*units::cm);
    ps.push_back(p);
    
    g1->SetPoint(i,x->at(i),y->at(i),z->at(i));
  }
  pcloud.AddPoints(ps);
  pcloud.build_kdtree_index();
  ps.clear();
  
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
  Double_t dQ1,dx1;
  T_rec->SetBranchAddress("x",&x1);
  T_rec->SetBranchAddress("y",&y1);
  T_rec->SetBranchAddress("z",&z1);
  T_rec->SetBranchAddress("q",&dQ1);
  T_rec->SetBranchAddress("nq",&dx1);
  TGraph *g2_xy = new TGraph();
  TGraph *g2_xz = new TGraph();
  TGraph *g2_yz = new TGraph();

  TGraph2D *g2 = new TGraph2D();

  TTree *t1 = new TTree("T","T");
  t1->SetDirectory(file);
  double max_dis=0; // maximum distance along the track 
  double beg_dis; // distance at the beginning ... 
  double end_dis; // distance at the end ...
  double total_dis2=0; // total distance^2
  double total_L = 0;
  int Npoints = 0;
  double total_dtheta = 0;
  double max_dtheta = 0;
  t1->Branch("max_dis",&max_dis);
  t1->Branch("beg_dis",&beg_dis);
  t1->Branch("end_dis",&end_dis);
  t1->Branch("total_dis2",&total_dis2);
  t1->Branch("total_L",&total_L);
  t1->Branch("N",&Npoints);
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
  
  t1->Branch("x",&x2);
  t1->Branch("y",&y2);
  t1->Branch("z",&z2);
  t1->Branch("dQ_rec",&dQ_rec);
  t1->Branch("dQ_tru",&dQ_tru);
  t1->Branch("dx",&dx);
  
  t1->Branch("x_pair",&x2_pair);
  t1->Branch("y_pair",&y2_pair);
  t1->Branch("z_pair",&z2_pair);

  t1->Branch("dis",&dis);
  t1->Branch("L",&L);
  t1->Branch("dtheta",&dtheta);

  double temp_L = 0;
  double prev_x1, prev_y1, prev_z1;

  Npoints = T_rec->GetEntries();

  std::map<std::tuple<int,int,int>, int> map_point_index;
  
  for (int i=0;i!=T_rec->GetEntries();i++){
    T_rec->GetEntry(i);
    
    
    if (i!=0){
      temp_L += sqrt(pow(x1-prev_x1,2)+pow(y1-prev_y1,2)+pow(z1-prev_z1,2));
    }
    
    // binning effect 1 us later from the binned slice effect, rebinned 4 ...
    // speed of imaging 1.101 mm / us
    // speed of simulation 1.098 mm/us
    // there is a potential 1 us offset at SP
    // -0.6 cm is the distance difference between Y and U planes
    x1 = (x1+0.1101)/1.1009999*1.098-0.6;//+ 4*0.1101; 
    
    g2_xy->SetPoint(i,x1,y1);
    g2_xz->SetPoint(i,x1,z1);
    g2_yz->SetPoint(i,y1,z1);

    Point p(x1*units::cm, y1*units::cm, z1*units::cm);
    ps.push_back(p);

    std::pair<double, Point> point_pair = pcloud.get_closest_point(p);

    map_point_index[std::make_tuple(p.x/(0.01*units::mm),p.y/(0.01*units::mm),p.z/(0.01*units::mm))] = x2->size();
    
    x2->push_back(x1);
    y2->push_back(y1);
    z2->push_back(z1);
    dQ_rec->push_back(dQ1);
    dQ_tru->push_back(0);
    dx->push_back(dx1);
    dis->push_back(point_pair.first/units::cm);

    x2_pair->push_back(point_pair.second.x/units::cm);
    y2_pair->push_back(point_pair.second.y/units::cm);
    z2_pair->push_back(point_pair.second.z/units::cm);
    
    L->push_back(temp_L);
    
    if (max_dis < point_pair.first/units::cm)
      max_dis = point_pair.first/units::cm;
    total_dis2 += pow(point_pair.first/units::cm,2);

    if (i==0){
      double dis1 = pow(x1-x->front(),2) + pow(y1-y->front(),2) + pow(z1-z->front(),2);
      double dis2 = pow(x1-x->back(),2) + pow(y1-y->back(),2) + pow(z1-z->back(),2);
      
      //std::cout << sqrt(dis1)/units::cm << " " << sqrt(dis2)/units::cm << std::endl;
      if (dis1 < dis2){
	beg_dis = sqrt(dis1);
      }else{
	beg_dis = sqrt(dis2);
      }
    }
    if (i==T_rec->GetEntries()-1){
      double dis1 = pow(x1-x->front(),2) + pow(y1-y->front(),2) + pow(z1-z->front(),2);
      double dis2 = pow(x1-x->back(),2) + pow(y1-y->back(),2) + pow(z1-z->back(),2);
      if (dis1 < dis2){
	end_dis = sqrt(dis1);
      }else{
	end_dis = sqrt(dis2);
      }
    }
    

    
    //  std::cout << point_pair.first/units::cm << " " << point_pair.second.x/units::cm << " " << point_pair.second.y/units::cm << " " << point_pair.second.z/units::cm << std::endl;
    
    g2->SetPoint(i,x1,y1,z1);

    prev_x1 = x1;
    prev_y1 = y1;
    prev_z1 = z1;
    
  }
  
  pcloud1.AddPoints(ps);
  pcloud1.build_kdtree_index();

   for (size_t i=0;i!=x->size();i++){
     Point p(x->at(i)*units::cm,y->at(i)*units::cm,z->at(i)*units::cm);
     std::pair<double, Point> point_pair = pcloud1.get_closest_point(p);
     int index = map_point_index[std::make_tuple(int(point_pair.second.x/(0.01*units::mm))
						 ,int(point_pair.second.y/(0.01*units::mm)),int(point_pair.second.z/(0.01*units::mm)))];
     //std::cout << index << std::endl;
     dQ_tru->at(index) += Q->at(i);

     if (i==0 || i==x->size()-1)
       std::cout << p << " " << sqrt(pow(p.x-point_pair.second.x,2)+pow(p.y-point_pair.second.y,2)+pow(p.z-point_pair.second.z,2))/units::cm << std::endl;
     
     // g1->SetPoint(i,x->at(i),y->at(i),z->at(i));
  }
  

  if (x2->size()>1){
    for (size_t i=0;i!=x2->size();i++){
      if (i==0){
	TVector3 dir1(x2->at(1)-x2->at(0),
		      y2->at(1)-y2->at(0),
		      z2->at(1)-z2->at(0));
	TVector3 dir2(x2_pair->at(1) - x2_pair->at(0),
		      y2_pair->at(1) - y2_pair->at(0),
		      z2_pair->at(1) - z2_pair->at(0));
	dtheta->push_back(dir1.Angle(dir2));
      }else if(i==x2->size()-1){
	TVector3 dir1(x2->at(x2->size()-1) - x2->at(x2->size()-2),
		      y2->at(x2->size()-1) - y2->at(x2->size()-2),
		      z2->at(x2->size()-1) - z2->at(x2->size()-2));
	TVector3 dir2(x2_pair->at(x2->size()-1) - x2_pair->at(x2->size()-2),
		      y2_pair->at(x2->size()-1) - y2_pair->at(x2->size()-2),
		      z2_pair->at(x2->size()-1) - z2_pair->at(x2->size()-2));
	dtheta->push_back(dir1.Angle(dir2));
      }else{
	TVector3 dir1(x2->at(i+1)-x2->at(i),
		      y2->at(i+1)-y2->at(i),
		      z2->at(i+1)-z2->at(i));
	TVector3 dir2(x2_pair->at(i+1) - x2_pair->at(i),
		      y2_pair->at(i+1) - y2_pair->at(i),
		      z2_pair->at(i+1) - z2_pair->at(i));

	TVector3 dir3(x2->at(i-1)-x2->at(i),
		      y2->at(i-1)-y2->at(i),
		      z2->at(i-1)-z2->at(i));
	TVector3 dir4(x2_pair->at(i-1) - x2_pair->at(i),
		      y2_pair->at(i-1) - y2_pair->at(i),
		      z2_pair->at(i-1) - z2_pair->at(i));
	dtheta->push_back((dir1.Angle(dir2)+dir3.Angle(dir4))/2.);
      }
      if (dtheta->back() > max_dtheta)
	max_dtheta = dtheta->back();
      total_dtheta += dtheta->back();
    }
  }else{
    dtheta->push_back(0);
  }
  
  total_L = temp_L;
  
  t1->Fill();
  
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
  t1->Write();
  

  // T_true->CloneTree(-1,"fast");
  // T_rec->CloneTree(-1,"fast");
  
  file->Write();
  file->Close();
  
  
}
