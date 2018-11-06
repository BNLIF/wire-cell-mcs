// Created on 10/31/2018
// by Wenqiang Gu (wgu@bnl.gov)
// to compile it standalone: g++ ToyMCS.cxx -o ToyMCS.exe `root-config --cflags --libs`

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

#include "WireCellMCSSim/MCSSim.h"


using namespace WireCell;

// c.f. http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf


int main(){
  MCStrack atrack;

  // particle type ...
  // 1: pi, 2 kaon, 3 muon, 4 proton 5, electron
  int particle_type = 3;

  double T_init = 100*units::MeV; // initial kinetic energy
  TVector3 pos_init(1.5*units::m,30*units::cm,6*units::m);
  TVector3 dir_init(0.4,0.4,0.4);
  dir_init *= 1./dir_init.Mag();

  double step_size = 0.1*units::cm;

  MCSTrackSim(atrack, particle_type, T_init, pos_init, dir_init, step_size);

  TFile* ofile = new TFile("mcs-tracks.root","RECREATE");
  TTree* T = new TTree("T","tracks and vertices");
  T->Branch("N", &atrack.N, "N/I");
  T->Branch("x", &atrack.x);
  T->Branch("y", &atrack.y);
  T->Branch("z", &atrack.z);
  T->Branch("Q", &atrack.Q);
  T->SetDirectory(ofile);
  
  // make out file ...
  std::ofstream outfile("mcssim.json");
  outfile << "{" << std::endl;
  outfile << " \"depos\": [" << std::endl;

  for (int i=0;i!=atrack.N;i++){
    outfile << "  {" << std::endl;

    outfile << "   \"n\": " << atrack.Q.at(i) << "," << std::endl;
    outfile << "   \"q\": 0," << std::endl;
    outfile << "   \"s\": 0," << std::endl;
    outfile << "   \"t\": 0," << std::endl;
    outfile << "   \"y\": " << atrack.y.at(i)/units::mm << "," << std::endl;
    outfile << "   \"x\": " << atrack.x.at(i)/units::mm << "," << std::endl;
    outfile << "   \"z\": " << atrack.z.at(i)/units::mm  << std::endl;

    outfile << "  }," << std::endl;
    
    atrack.y.at(i)/=units::cm;
    atrack.x.at(i)/=units::cm;
    atrack.z.at(i)/=units::cm;
    
    
  }

  
  outfile << "  {" << std::endl;
  outfile << "  }" << std::endl;
  outfile << " ]" << std::endl;
  outfile << "}" << std::endl;
  outfile.close();

  T->Fill();
  T->Write();
  ofile->Close();
  
  // TFile* ofile = new TFile("mcs-tracks.root","RECREATE");
  // TTree* T = new TTree("T","tracks and vertices");
  // T->Branch("N", &atrack.N, "N/I");
  // T->Branch("x", &atrack.x);
  // T->Branch("y", &atrack.y);
  // T->Branch("z", &atrack.z);

  // for(int i=0; i<1000/*ntracks*/; i++){
  //   double x = gRandom->Uniform(2.56*units::m);
  //   double y = gRandom->Uniform(2.3*units::m);
  //   double z = gRandom->Uniform(10.4*units::m);
  //   double cosTheta = gRandom->Uniform(-1,1);
  //   double phi = gRandom->Uniform(0,2*3.1415926);
  //   std::vector<double> initpos = {x,y,z};
  //   std::vector<double> initdir = {acos(cosTheta), phi};

  //   MultiScattSim(atrack, 100/*nvertices*/, initpos, initdir);

  //   for (size_t j=0; j!=atrack.x.size(); j++){
  //     atrack.x.at(j) /= units::m;
  //     atrack.y.at(j) /= units::m;
  //     atrack.z.at(j) /= units::m;
  //   }
    
  //   T->Fill();
  // }
  // T->Write();
  // ofile->Close();
}
