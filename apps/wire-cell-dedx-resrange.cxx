#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"

#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

#include "WCPMCSSim/Eloss.h"

using namespace WCP;
using namespace std;

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}


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
  auto hdedx_mpv_resrange = new TH2F("hdEdxMPVRR","",1200,0,120, 600,0,30); 
  auto hdedx_T = new TH2F("hdedxT","",1000,0,10,250,2,4.5);

  std::vector<double> mean_stoppowers, mpv_stoppowers, resranges; // save the dEdx and residual range

  for(int iter=0; iter<num_particles; iter++) {
    if (iter%10==0) {
      std::cout << "===" << std::endl;
      std::cout << "=" << std::endl;
      std::cout << "= iteration of particles: " << iter << std::endl;
      std::cout << "=" << std::endl;
      std::cout << "===" << std::endl;
    }
    double T = init_ke;
    std::vector<double> edeps, Tcollection;
    double  last_T = 0;
    while(T>0){
      // double dEdx_MPV = loss.get_MPV_dEdx(T, step_size);
      
      // double dEdx_mean_ioniz = loss.get_mean_ioniz_dEdx(T);
      // double dEdx_mean_total = loss.get_mean_total_dEdx(T);
  
      double dEdx_mean_dEdx = loss.get_mean_dEdx(T);
      // std::cout << "T[MeV]: " << T/units::MeV << " dEdx_mean: " << dEdx_mean_dEdx << std::endl;
      // hdedx_T->Fill(T/units::GeV, dEdx_mean_dEdx/(units::MeV/units::cm));

      double eloss = step_size*loss.get_dEdx(T, step_size);
      // double eloss = step_size*dEdx_mean_dEdx;


      if (T>eloss) {
        edeps.push_back(eloss);
        Tcollection.push_back(T);
        T -= eloss;
      }
      else {
        last_T = T;
        T = 0;
      }
    }

    // special treatment for the last step
    // double dEdx1 = edeps.back()/step_size;
    // double dEdx2 = loss.get_MPV_dEdx(last_T, step_size); // mean or MPV?
    // double dEdx_avg = 0.5*(dEdx1 + dEdx2);
    double dEdx_avg = loss.get_mean_dEdx(last_T);
    double last_step_length = last_T/dEdx_avg;
    std::cout << "T: " << last_T/units::MeV << " last step length: " << last_step_length/units::cm << std::endl;


    int istep=0;
    double total_length = edeps.size() * step_size + last_step_length;
    for(auto dE: edeps) {
      double dEdx = dE/step_size;
      double dEdx_mpv = loss.get_MPV_dEdx(Tcollection.at(istep), step_size);
      double dEdx_mean = loss.get_mean_dEdx(Tcollection.at(istep)); // get_mean_total_dEdx if input file is given
      istep ++;
      double res_range = total_length - step_size*istep;
      hdedx_resrange->Fill(res_range/units::cm, dEdx/(units::MeV/units::cm));
      hdedx_mpv_resrange->Fill(res_range/units::cm, dEdx_mpv/(units::MeV/units::cm));

      mpv_stoppowers.push_back(dEdx_mpv/(units::MeV/units::cm));
      mean_stoppowers.push_back(dEdx_mean/(units::MeV/units::cm));
      resranges.push_back(res_range/units::cm);
    }
  }

  ofstream ofs("mp_dedx_mean_vs_range_1GeV_WCP.txt", std::ofstream::out);
  for (auto ind: sort_indexes(resranges)) {
    ofs << resranges[ind] << " " << mean_stoppowers[ind] << endl;
  }
  ofs.close();

  ofstream ofs1("mp_dedx_mpv_vs_range_1GeV_WCP.txt", std::ofstream::out);
  for (auto ind: sort_indexes(resranges)) {
    ofs1 << resranges[ind] << " " << mpv_stoppowers[ind] << endl;
  }
  ofs1.close();

  auto ofile = new TFile("ofile.root","recreate");
  hdedx_resrange->Write();
  hdedx_mpv_resrange->Write();
  hdedx_T->Write();
  ofile->Close();

}
