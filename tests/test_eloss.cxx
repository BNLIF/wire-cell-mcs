#include "WireCellMCSSim/Eloss.h"
#include "WireCellData/LAr.h"

#include <iostream>

using namespace WireCell;
using namespace std;




int main(int argc, char* argv[])
{
  Eloss cal_eloss(4); // pion, kaon, muon, proton

  std::cout << cal_eloss.get_mean_dEdx(100*units::MeV) / units::MeV * units::cm << std::endl;

  LAr lar;
  //  std::cout << lar.Ldensity(89) << std::endl;

  lar.print_critical();
  lar.print_triple();
  lar.print_boiling();

  std::cout << lar.Gdensity(89*units::kelvin)/(units::g/pow(units::cm,3)) << " " << lar.Ldensity(89*units::kelvin)/(units::g/pow(units::cm,3)) << " " << lar.VPressure(89*units::kelvin)/units::bar << " " << lar.vDarIon(89*units::kelvin,0.5*units::kilovolt/units::cm)/(units::m/units::second) << " " << lar.vDarIonV(89*units::kelvin,0.5*units::kilovolt/units::cm)/(units::m/units::second) <<std::endl;
  std::cout << lar.Viscosity(89*units::kelvin)/(1e-6 * units::pascal * units::second) << " " << lar.SpeedofSound(89*units::kelvin)/(units::m/units::second) << " " << lar.Cp(89*units::kelvin)/(1000 * units::joule/ units::kilogram /units::kelvin) << " " << lar.Cv(89*units::kelvin)/(1000 * units::joule/ units::kilogram /units::kelvin) << " " << lar.IsoCom(89*units::kelvin)/(pow(units::cm,2)/units::g/units::cm*pow(units::second,2)) << std::endl;
  std::cout << lar.Enthalpy(89*units::kelvin)/(1000*units::joule/units::mole) << " " << lar.epsilon(89*units::kelvin) << " " << lar.GInrf(150*units::nm) << " " << lar.LInrf(150*units::nm,89*units::kelvin) << " " << lar.RRLAr(150*units::nm,89*units::kelvin)/units::cm << std::endl;
  std::cout << lar.vD(89*units::kelvin,0.5*units::kilovolt/units::cm,1)/(units::mm/units::microsecond) << " " << lar.vDrift(89*units::kelvin,0.1*units::kilovolt/units::cm)/(units::mm/units::microsecond)<< " " << lar.Diffusion(89*units::kelvin, 0.5*units::kilovolt/units::cm,2)/units::eV << " " << lar.recombine_Birks(1.8*units::MeV/units::cm*lar.Ldensity(89*units::kelvin),89*units::kelvin,0.5*units::kilovolt/units::cm) << " " << lar.recombine_Box(1.8*units::MeV/units::cm*lar.Ldensity(89*units::kelvin),89*units::kelvin,0.5*units::kilovolt/units::cm) << std::endl;
  std::cout << lar.SubLimeAr(82*units::kelvin)/units::bar << " " << lar.MeltAr(82*units::kelvin)/units::bar << " " << lar.BoilAr(82*units::kelvin)/units::bar << " " << lar.ele_lifetime(1,89*units::kelvin, 0.5*units::kilovolt/units::cm, 5)/units::millisecond << std::endl;
  return 0;
  //  std::cout << "test" << std::endl;
}
