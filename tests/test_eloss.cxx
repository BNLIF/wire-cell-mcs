#include "WireCellMCSSim/Eloss.h"

#include <iostream>

using namespace WireCell;
using namespace std;




int main(int argc, char* argv[])
{
  Eloss cal_eloss(4); // pion, kaon, muon, proton

  std::cout << cal_eloss.get_mean_dEdx(100*units::MeV) / units::MeV * units::cm << std::endl;
  
  return 0;
  //  std::cout << "test" << std::endl;
}
