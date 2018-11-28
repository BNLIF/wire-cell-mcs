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

using namespace MCS;

// c.f. http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf


int main(int argc, char* argv[]){
  if (argc < 2){
    std::cerr << "wire-cell-uboone -p[particle type] -t[kinetic energy] -x[x-pos] -y[y-pos] -z[z-pos] -a[x-dir] -b[y dir] -c[z dir] -s[step size]" << std::endl;
    std::cerr << "particle type:  1 pi, 2 kaon, 3 muon, 4 proton, 5 electron" << std::endl;
    std::cerr << "kinetic energy in MeV" << std::endl;
    std::cerr << "position in cm" << std::endl;
    std::cerr << "step size in cm" << std::endl;
    
    return 1;
  }
  
  

  // particle type ...
  // 1: pi, 2 kaon, 3 muon, 4 proton 5, electron
  int particle_type = 5;
  double T_init = 0.5*units::MeV; // initial kinetic energy
  TVector3 pos_init(1.5*units::m,30*units::cm,6*units::m);
  TVector3 dir_init(0.4,0.4,0.4);
  double step_size = 0.1*units::cm;
  
  for(Int_t i = 1; i != argc; i++){
    switch(argv[i][1]){
    case 'p':
      particle_type = atoi(&argv[i][2]); 
      break;
    case 't':
      T_init = atof(&argv[i][2])*units::MeV; 
      break;
    case 'x':
      pos_init.SetX(atof(&argv[i][2])*units::cm); 
      break;
    case 'y':
      pos_init.SetY(atof(&argv[i][2])*units::cm);
      break;
    case 'z':
      pos_init.SetZ(atof(&argv[i][2])*units::cm);
      break;
    case 'a':
      dir_init.SetX(atof(&argv[i][2]));
      break;
    case 'b':
      dir_init.SetY(atof(&argv[i][2]));
      break;
    case 'c':
      dir_init.SetZ(atof(&argv[i][2]));
      break;
    case 's':
      step_size = atof(&argv[i][2])*units::cm;
    }
  }
  

  
  dir_init *= 1./dir_init.Mag();

  MCStrack atrack;
  
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

  for (int j=0;j!=50;j++){
    TVector3 pos_gen;
    pos_gen.SetXYZ(pos_init.X()+gRandom->Gaus(0,25*units::cm),
		   pos_init.Y()+gRandom->Gaus(0,25*units::cm),
		   pos_init.Z()+gRandom->Gaus(0,25*units::cm));
    
    
    MCSTrackSim(atrack, particle_type, T_init, pos_gen, dir_init, step_size);
    
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
    
    T->Fill();
    atrack.clear();
  }
  
  outfile << "  {" << std::endl;
  outfile << "  }" << std::endl;
  outfile << " ]" << std::endl;
  outfile << "}" << std::endl;
  outfile.close();

  
  T->Write();
  ofile->Close();
  
 
}
