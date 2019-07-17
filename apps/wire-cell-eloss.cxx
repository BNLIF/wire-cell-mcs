#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

#include "WireCellMCSSim/MCSSim.h"
#include "WireCellMCSSim/Eloss.h"

using namespace WireCell;

using namespace MCS;

// c.f. http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf


int main(int argc, char* argv[]){
  int particle_type = 2; //kaon ...
  Eloss loss(2);
  double T = 10000*units::MeV;
  double L = 0;
  int counter = 0;
  
  while(T>0){
    double dEdx_MPV = loss.get_MPV_dEdx(T, 3*units::mm);
    double dEdx_mean_ioniz = loss.get_mean_ioniz_dEdx(T);
    double dEdx_mean_total = loss.get_mean_total_dEdx(T);

    L += 0.3*units::mm;
    T -= 0.3*units::mm * dEdx_mean_total;
    if (counter%5==0)
      std::cout << 387.33 + 3844.35 - L/units::cm << " " << T/units::MeV << " " << dEdx_MPV/(units::MeV/units::cm) << " " << dEdx_mean_ioniz/(units::MeV/units::cm) << " " << dEdx_mean_total/(units::MeV/units::cm) << std::endl;
  }
  
}
