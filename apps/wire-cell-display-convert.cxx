
#include "WCPData/ToyPointCloud.h"

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

using namespace WCP;

int main(int argc, char* argv[])
{
  if (argc < 2){
    std::cerr << "usage: wire-cell-track-com -a[truth.root] -b[reco.root] -t[reco_treename] -n[truth_treename] -o[out.root] -f[1:MC,2:data]" << std::endl;
    return 1;
  }
  
  TString reco_filename = "tracking_0_0_0.root"; 
  TString reco_treename = "T_rec_charge";
  TString proj_treename = "T_proj_data";
  
  TString truth_filename = "mcs-tracks.root"; 
  TString truth_treename = "T";
  TString out_filename = "track_com.root";

  int file_type = 1; // 1 for MC and 2 for data ...

 
  
  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 'b':
       reco_filename = &argv[i][2];
       break;
     case 't':
       reco_treename = &argv[i][2];
       break;
     case 'a':
       truth_filename = &argv[i][2];
       break;
     case 'n':
       truth_treename = &argv[i][2]; 
       break;
     case 'o':
       out_filename = &argv[i][2];
       break;
     case 'f':
       file_type = atoi(&argv[i][2]);
       break;
     }
  }

  if (file_type==1){
    std::cout << truth_filename << " " << reco_filename << " " << truth_treename << " " << reco_treename << std::endl;
  }else{
    std::cout << reco_filename << " " << reco_treename << std::endl;
  }

  TFile *file1 = new TFile(reco_filename);
  TTree *T_bad_ch = (TTree*)file1->Get("T_bad_ch");
   
  
  TFile *file = new TFile(out_filename,"RECREATE");
  if (T_bad_ch!=0){
    T_bad_ch->CloneTree()->Write();
  }

  std::vector<std::vector<double> > *vx = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *vy = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *vz = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *vQ = new std::vector<std::vector<double> >;  
  std::vector<int> *vN = new std::vector<int>;

  std::vector<double> *x = new std::vector<double>;
  std::vector<double> *y = new std::vector<double>;
  std::vector<double> *z = new std::vector<double>;
  std::vector<double> *Q = new std::vector<double>;
  Int_t N;
  
  if (file_type==1){ // MC data ... 
    TChain *T_true = new TChain(truth_treename,truth_treename);
    T_true->Add(truth_filename);
    TTree *t2 = new TTree("T_true","T_true");
    T_true->SetBranchAddress("N",&N);
    T_true->SetBranchAddress("x",&x);
    T_true->SetBranchAddress("y",&y);
    T_true->SetBranchAddress("z",&z);
    T_true->SetBranchAddress("Q",&Q);
    T_true->GetEntry(0);
    vN->push_back(N);
    vx->push_back(*x);
    vy->push_back(*y);
    vz->push_back(*z);
    vQ->push_back(*Q);

    for (size_t i=0;i!=vx->at(0).size();i++){
      vx->at(0).at(i) = (vx->at(0).at(i)+0.6)/1.098*1.1009999-0.1101;
    }
    
    t2->Branch("N",&vN);
    t2->Branch("x",&vx);
    t2->Branch("y",&vy);
    t2->Branch("z",&vz);
    t2->Branch("Q",&vQ);
    t2->Fill();
  }

  
  
  // std::cout << N << std::endl;

  // TGraph *g1_xy = new TGraph();
  // TGraph *g1_xz = new TGraph();
  // TGraph *g1_yz = new TGraph();

  // TGraph2D *g1 = new TGraph2D();

  ToyPointCloud pcloud, pcloud1;

  WCP::PointVector ps;
  
  if (file_type==1){ // MC ... 
    for (size_t i=0;i!=x->size();i++){
      // g1_xy->SetPoint(i,x->at(i),y->at(i));
      // g1_xz->SetPoint(i,x->at(i),z->at(i));
      // g1_yz->SetPoint(i,y->at(i),z->at(i));

      x->at(i) = (x->at(i)+0.6)/1.098*1.1009999-0.1101;
      //      x1 = (x1+0.1101)/1.1009999*1.098-0.6;//+ 4*0.1101; 
      
      Point p(x->at(i)*units::cm,y->at(i)*units::cm,z->at(i)*units::cm);
      ps.push_back(p);
      //g1->SetPoint(i,x->at(i),y->at(i),z->at(i));
    }
    pcloud.AddPoints(ps);
    pcloud.build_kdtree_index();
  }
  
  ps.clear();
  
  
  // g1_xy->SetLineColor(1);
  // g1_xz->SetLineColor(1);
  // g1_yz->SetLineColor(1);

  // g1_xy->SetLineWidth(2);
  // g1_xz->SetLineWidth(2);
  // g1_yz->SetLineWidth(2);

  // g1->SetLineColor(1);
  // g1->SetLineWidth(2);

 
  
  TChain *T_rec = new TChain(reco_treename,reco_treename);
  TChain *T_proj_data = new TChain(proj_treename,proj_treename);
  TChain *T_proj = new TChain("T_proj","T_proj");
  
  T_rec->Add(reco_filename);
  T_proj_data->Add(reco_filename);
  T_proj->Add(reco_filename);
  Double_t x1,y1,z1;
  Double_t dQ1,dx1,ndf;
  Double_t pu, pv, pw, pt;

  T_rec->SetBranchAddress("x",&x1);
  T_rec->SetBranchAddress("y",&y1);
  T_rec->SetBranchAddress("z",&z1);
  T_rec->SetBranchAddress("q",&dQ1);
  T_rec->SetBranchAddress("nq",&dx1);
  T_rec->SetBranchAddress("ndf",&ndf);
  T_rec->SetBranchAddress("pu",&pu);
  T_rec->SetBranchAddress("pv",&pv);
  T_rec->SetBranchAddress("pw",&pw);
  T_rec->SetBranchAddress("pt",&pt);

  Double_t reduced_chi2;
  if (T_rec->GetBranch("reduced_chi2")){
    T_rec->SetBranchAddress("reduced_chi2",&reduced_chi2);
  }
  
  TTree *t1 = new TTree("T_rec","T_rec");
  t1->SetDirectory(file);
  std::vector<double> *max_dis = new std::vector<double>; // maximum distance along the track 
  std::vector<double> *beg_dis = new std::vector<double>; // distance at the beginning ... 
  std::vector<double> *end_dis = new std::vector<double>; // distance at the end ...
  std::vector<double> *total_dis2 = new std::vector<double>; // total distance^2
  std::vector<double> *total_L = new std::vector<double>;
  std::vector<int> *Npoints = new std::vector<int>;
  std::vector<double> *total_dtheta = new std::vector<double>;
  std::vector<double> *max_dtheta = new std::vector<double>;
  if (file_type==1){
    t1->Branch("stat_max_dis",&max_dis);
    t1->Branch("stat_beg_dis",&beg_dis);
    t1->Branch("stat_end_dis",&end_dis);
    t1->Branch("stat_total_dis2",&total_dis2);
    t1->Branch("stat_total_dtheta",&total_dtheta);
    t1->Branch("stat_max_dtheta",&max_dtheta);
  }
  t1->Branch("stat_total_L",&total_L);
  t1->Branch("stat_N",&Npoints);
  
  
  std::vector<std::vector<double> > *x2 = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *y2 = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *z2 = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *dQ_rec = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *dQ_tru = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *dx = new std::vector<std::vector<double> >;
  std::vector<int> *cluster_id = new std::vector<int>;
  std::vector<std::vector<double> > *rec_pu = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *rec_pv = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *rec_pw = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *rec_pt = new std::vector<std::vector<double> >;
  
  std::vector<std::vector<double> > *x2_pair = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *y2_pair = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *z2_pair = new std::vector<std::vector<double> >;

  std::vector<std::vector<double> > *dis = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *L = new std::vector<std::vector<double> >;
  std::vector<std::vector<double> > *dtheta = new std::vector<std::vector<double> >;

  std::vector<std::vector<double> > *Vreduced_chi2 = new std::vector<std::vector<double> > ;

  
  t1->Branch("rec_x",&x2);
  t1->Branch("rec_y",&y2);
  t1->Branch("rec_z",&z2);
  t1->Branch("rec_dQ",&dQ_rec);
  t1->Branch("rec_dx",&dx);
  t1->Branch("rec_L",&L);
  t1->Branch("rec_cluster_id",&cluster_id);
  t1->Branch("rec_u",&rec_pu);
  t1->Branch("rec_v",&rec_pv);
  t1->Branch("rec_w",&rec_pw);
  t1->Branch("rec_t",&rec_pt);
  if (T_rec->GetBranch("reduced_chi2")){
    t1->Branch("reduced_chi2",&Vreduced_chi2);
  }
  
  
  if (file_type==1){
    t1->Branch("true_dQ",&dQ_tru);
    t1->Branch("true_x",&x2_pair);
    t1->Branch("true_y",&y2_pair);
    t1->Branch("true_z",&z2_pair);
    
    t1->Branch("com_dis",&dis);
    t1->Branch("com_dtheta",&dtheta);
  }
  
  double temp_L = 0;
  double prev_x1, prev_y1, prev_z1;
  int prev_cluster_id = -1;

  // Npoints = T_rec->GetEntries();

  std::map<std::tuple<int,int,int>, std::pair<int, int> > map_point_index;

  //  std::cout << T_rec->GetEntries() << std::endl;
  
  for (int i=0;i!=T_rec->GetEntries();i++){
    T_rec->GetEntry(i);
    if (std::round(ndf)!=prev_cluster_id){
      Npoints->push_back(0);
      total_L->push_back(0);

      std::vector<double> temp_dQ_tru;
      dQ_tru->push_back(temp_dQ_tru);
      std::vector<double> temp_dis;
      dis->push_back(temp_dis);
      std::vector<double> temp_x2_pair;
      x2_pair->push_back(temp_x2_pair);
      std::vector<double> temp_y2_pair;
      y2_pair->push_back(temp_y2_pair);
      std::vector<double> temp_z2_pair;
      z2_pair->push_back(temp_z2_pair);

      std::vector<double> temp_x2;
      x2->push_back(temp_x2);
      std::vector<double> temp_y2;
      y2->push_back(temp_y2);
      std::vector<double> temp_z2;
      z2->push_back(temp_z2);
      std::vector<double> temp_rec_pu;
      rec_pu->push_back(temp_rec_pu);
      std::vector<double> temp_rec_pv;
      rec_pv->push_back(temp_rec_pv);
      std::vector<double> temp_rec_pw;
      rec_pw->push_back(temp_rec_pw);
      std::vector<double> temp_rec_pt;
      rec_pt->push_back(temp_rec_pt);
      cluster_id->push_back(std::round(ndf));
      std::vector<double> temp_dQ_rec;
      dQ_rec->push_back(temp_dQ_rec);
      std::vector<double> temp_dx;
      dx->push_back(temp_dx);
      std::vector<double> temp_L;
      L->push_back(temp_L);

      if (T_rec->GetBranch("reduced_chi2")){
	std::vector<double> temp_reduced_chi2;
	Vreduced_chi2->push_back(temp_reduced_chi2);
      }
      
      max_dis->push_back(0);
      total_dis2->push_back(0);
    }
    Npoints->back()++;
    
    if (Npoints->back()!=1){
      total_L->back() += sqrt(pow(x1-prev_x1,2)+pow(y1-prev_y1,2)+pow(z1-prev_z1,2));
    }
    
    // binning effect 1 us later from the binned slice effect, rebinned 4 ...
    // speed of imaging 1.101 mm / us
    // speed of simulation 1.098 mm/us
    // there is a potential 1 us offset at SP
    // -0.6 cm is the distance difference between Y and U planes
    //    x1 = (x1+0.1101)/1.1009999*1.098-0.6;//+ 4*0.1101; 
    Point p(x1*units::cm, y1*units::cm, z1*units::cm);
    ps.push_back(p);

   
    
    if (file_type==1){
      std::pair<double, Point> point_pair = pcloud.get_closest_point(p);
      map_point_index[std::make_tuple(p.x/(0.01*units::mm),p.y/(0.01*units::mm),p.z/(0.01*units::mm))] = std::make_pair(x2->size()-1,x2->back().size());
      dQ_tru->back().push_back(0);
      dis->back().push_back(point_pair.first/units::cm);
      

      x2_pair->back().push_back(point_pair.second.x/units::cm);
      y2_pair->back().push_back(point_pair.second.y/units::cm);
      z2_pair->back().push_back(point_pair.second.z/units::cm);
        
      if (max_dis->back() < point_pair.first/units::cm)
     	max_dis->back() = point_pair.first/units::cm;
      total_dis2->back() += pow(point_pair.first/units::cm,2);

      if (Npoints->back()==1){
     	double dis1 = pow(x1-x->front(),2) + pow(y1-y->front(),2) + pow(z1-z->front(),2);
     	double dis2 = pow(x1-x->back(),2) + pow(y1-y->back(),2) + pow(z1-z->back(),2);
      
    // 	//std::cout << sqrt(dis1)/units::cm << " " << sqrt(dis2)/units::cm << std::endl;
     	if (dis1 < dis2){
     	  beg_dis->push_back(sqrt(dis1));
     	}else{
	  beg_dis->push_back(sqrt(dis2));
     	}
	end_dis->push_back(0);
      }
    //   if (i==T_rec->GetEntries()-1)
      {
     	double dis1 = pow(x1-x->front(),2) + pow(y1-y->front(),2) + pow(z1-z->front(),2);
     	double dis2 = pow(x1-x->back(),2) + pow(y1-y->back(),2) + pow(z1-z->back(),2);
     	if (dis1 < dis2){
     	  end_dis->back() = sqrt(dis1);
     	}else{
     	  end_dis->back() = sqrt(dis2);
     	}
      }
    }

    x2->back().push_back(x1);
    y2->back().push_back(y1);
    z2->back().push_back(z1);
    rec_pu->back().push_back(pu);
    rec_pv->back().push_back(pv);
    rec_pw->back().push_back(pw);
    rec_pt->back().push_back(pt);
    dQ_rec->back().push_back(dQ1);
    dx->back().push_back(dx1);
    L->back().push_back(total_L->back());
    
    if (T_rec->GetBranch("reduced_chi2")){
      Vreduced_chi2->back().push_back(reduced_chi2);
    }
  

    prev_x1 = x1;
    prev_y1 = y1;
    prev_z1 = z1;
    prev_cluster_id = std::round(ndf);

    //std::cout << prev_cluster_id << std::endl;
  } // loop over i ...
  
  
  pcloud1.AddPoints(ps);
  pcloud1.build_kdtree_index();




  //std::cout << "haha " << std::endl;
  if (file_type==1){
    for (size_t i=0;i!=x->size();i++){
      Point p(x->at(i)*units::cm,y->at(i)*units::cm,z->at(i)*units::cm);
      std::pair<double, Point> point_pair = pcloud1.get_closest_point(p);
      int index = map_point_index[std::make_tuple(int(point_pair.second.x/(0.01*units::mm))
  						  ,int(point_pair.second.y/(0.01*units::mm)),int(point_pair.second.z/(0.01*units::mm)))].second;
      int index1 = map_point_index[std::make_tuple(int(point_pair.second.x/(0.01*units::mm))
  						  ,int(point_pair.second.y/(0.01*units::mm)),int(point_pair.second.z/(0.01*units::mm)))].first;

      // std::cout << index << " " << dQ_tru->back().size() << " " << dQ_tru->front().size() << std::endl;
      if (index <dQ_tru->at(index1).size())
	dQ_tru->at(index1).at(index) += Q->at(i);
      
      // if (i==0 || i==x->size()-1)
      //   std::cout << p << " " << sqrt(pow(p.x-point_pair.second.x,2)+pow(p.y-point_pair.second.y,2)+pow(p.z-point_pair.second.z,2))/units::cm << std::endl;
      // g1->SetPoint(i,x->at(i),y->at(i),z->at(i));
    }

    //    std::cout << x2->back().size() << std::endl;

    for (size_t k=0;k!=x2->size();k++){
      if (x2->at(k).size()>1){
	std::vector<double> temp_dtheta;
	dtheta->push_back(temp_dtheta);
	max_dtheta->push_back(0);
	total_dtheta->push_back(0);
	
	for (size_t i=0;i!=x2->at(k).size();i++){
	  if (i==0){
	    TVector3 dir1(x2->at(k).at(1)-x2->at(k).at(0),
			  y2->at(k).at(1)-y2->at(k).at(0),
			  z2->at(k).at(1)-z2->at(k).at(0));
	    TVector3 dir2(x2_pair->at(k).at(1) - x2_pair->at(k).at(0),
			  y2_pair->at(k).at(1) - y2_pair->at(k).at(0),
			  z2_pair->at(k).at(1) - z2_pair->at(k).at(0));
	    dtheta->at(k).push_back(dir1.Angle(dir2));
	  }else if(i==x2->at(k).size()-1){
	    TVector3 dir1(x2->at(k).at(x2->at(k).size()-1) - x2->at(k).at(x2->at(k).size()-2),
			  y2->at(k).at(x2->at(k).size()-1) - y2->at(k).at(x2->at(k).size()-2),
			  z2->at(k).at(x2->at(k).size()-1) - z2->at(k).at(x2->at(k).size()-2));
	    TVector3 dir2(x2_pair->at(k).at(x2->at(k).size()-1) - x2_pair->at(k).at(x2->at(k).size()-2),
			  y2_pair->at(k).at(x2->at(k).size()-1) - y2_pair->at(k).at(x2->at(k).size()-2),
			  z2_pair->at(k).at(x2->at(k).size()-1) - z2_pair->at(k).at(x2->at(k).size()-2));
	    dtheta->at(k).push_back(dir1.Angle(dir2));
	  }else{
	    TVector3 dir1(x2->at(k).at(i+1)-x2->at(k).at(i),
			  y2->at(k).at(i+1)-y2->at(k).at(i),
			  z2->at(k).at(i+1)-z2->at(k).at(i));
	    TVector3 dir2(x2_pair->at(k).at(i+1) - x2_pair->at(k).at(i),
			  y2_pair->at(k).at(i+1) - y2_pair->at(k).at(i),
			  z2_pair->at(k).at(i+1) - z2_pair->at(k).at(i));
	    
	    TVector3 dir3(x2->at(k).at(i-1)-x2->at(k).at(i),
			  y2->at(k).at(i-1)-y2->at(k).at(i),
			  z2->at(k).at(i-1)-z2->at(k).at(i));
	    TVector3 dir4(x2_pair->at(k).at(i-1) - x2_pair->at(k).at(i),
			  y2_pair->at(k).at(i-1) - y2_pair->at(k).at(i),
			  z2_pair->at(k).at(i-1) - z2_pair->at(k).at(i));
	    dtheta->at(k).push_back((dir1.Angle(dir2)+dir3.Angle(dir4))/2.);
	  }
	  if (dtheta->at(k).back() > max_dtheta->at(k))
	    max_dtheta->at(k) = dtheta->at(k).back();
	  total_dtheta->at(k) += dtheta->at(k).back();
	}
      }else{
	dtheta->at(k).push_back(0);
      }
    }
  }
  t1->Fill();
  
  T_proj_data->CloneTree(-1,"fast");
  T_proj->CloneTree(-1,"fast");
  
  file->Write();
  file->Close();
  
  
}
