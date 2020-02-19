#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"

#include "WCPMCSSim/Eloss.h"

using namespace WCP;

// c.f. http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf


int main(int argc, char* argv[]){

  if (argc < 2){
    std::cerr << "usage: wire-cell-dedx-resrange -p[particle_type: 1-5] -k[initial KE in GeV] -n[num_particles] -s[step size in cm]" << std::endl;
    std::cerr << "particle type: 1->pion, 2->kaon, 3->muon, 4->proton, 5->electron" << std::endl;
    return 1;
  }

  int particle_type = 4; //proton
  int num_particles = 10000;
  double step_size = 0.5*units::cm;
  double init_ke = 1.0*units::GeV;
 
  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 'p':
       particle_type = atoi(&argv[i][2]);
       break;
     case 'k':
       init_ke = atof(&argv[i][2])*units::GeV;
       break;
     case 'n':
       num_particles = atoi(&argv[i][2]);
       break;
     case 's':
       step_size = atof(&argv[i][2])*units::cm; 
       break;
     }
  }

  Eloss loss(particle_type);

  std::cout << "particle_type: " << particle_type << std::endl;
  std::cout << "initial KE [GeV]: " << init_ke/units::GeV << std::endl;
  std::cout << "num_particles: " << num_particles << std::endl;
  std::cout << "step size [cm]: " << step_size/units::cm << std::endl;


  auto hdedx_resrange = new TH2F("hdEdxRR","",240,0,120, 600,0,30);
  auto hdedx_T = new TH2F("hdedxT","",1000,0,10,250,2,4.5);

  for(int iter=0; iter<num_particles; iter++) {
    if (iter%10==0) {
      std::cout << "===" << std::endl;
      std::cout << "=" << std::endl;
      std::cout << "= iteration of particles: " << iter << std::endl;
      std::cout << "=" << std::endl;
      std::cout << "===" << std::endl;
    }
    double T = init_ke;
    std::vector<double> edeps;
    while(T>0){
      // double dEdx_MPV = loss.get_MPV_dEdx(T, step_size);
      
      // double dEdx_mean_ioniz = loss.get_mean_ioniz_dEdx(T);
      // double dEdx_mean_total = loss.get_mean_total_dEdx(T);
  
      double dEdx_mean_dEdx = loss.get_mean_dEdx(T);
      // std::cout << "T[MeV]: " << T/units::MeV << " dEdx_mean: " << dEdx_mean_dEdx << std::endl;
      hdedx_T->Fill(T/units::GeV, dEdx_mean_dEdx/(units::MeV/units::cm));

      double eloss = step_size*loss.get_dEdx(T, step_size);

      T -= eloss;
      edeps.push_back(eloss);
    }

    int istep=0;
    double total_length = edeps.size() * step_size;
    for(auto dE: edeps) {
      istep ++;
      double dEdx = dE/step_size;
      double res_range = total_length - step_size*istep;
      hdedx_resrange->Fill(res_range/units::cm, dEdx/(units::MeV/units::cm));
    }
  }

  auto ofile = new TFile("ofile.root","recreate");
  hdedx_resrange->Write();
  hdedx_T->Write();
  ofile->Close();

}
