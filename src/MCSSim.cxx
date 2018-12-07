#include "WireCellMCSSim/MCSSim.h"
#include "WireCellMCSSim/Eloss.h"
#include "WireCellData/LAr.h"

#include <iostream>

using namespace WireCell;

using namespace MCS;

void WireCell::MCS::RotateUz(TVector3& direction, TVector3& v1){
  TVector3 unit = direction.Unit();
  v1.RotateUz(unit); 
}

void WireCell::MCS::LineToyTrackSim(MCStrack& atrack, TVector3 pos_init, TVector3 dir_init, int nstep, double step_size, int charge){

  dir_init = dir_init.Unit();
  
  atrack.N = 0;
  for (int i=0;i!=nstep;i++){
    atrack.x.push_back(pos_init.X());
    atrack.y.push_back(pos_init.Y());
    atrack.z.push_back(pos_init.Z());

    atrack.vx.push_back(dir_init.X());
    atrack.vy.push_back(dir_init.Y());
    atrack.vz.push_back(dir_init.Z());

    atrack.p.push_back(charge);

    atrack.Q.push_back(charge);
    pos_init.SetXYZ(pos_init.X() - dir_init.X() * step_size,
		    pos_init.Y() - dir_init.Y() * step_size,
		    pos_init.Z() - dir_init.Z() * step_size);
    atrack.N ++;
  }
}


void WireCell::MCS::LineTrackSim(MCStrack& atrack, int particle_type, double T_init, TVector3 pos_init, TVector3 dir_init, double step_size){
  Eloss cal_loss(particle_type); // only for other particles ...
  LAr lar;

  // proper normalize ... 
  atrack.N = 0;
  dir_init *= 1./dir_init.Mag();
  
  while (T_init >0){
    atrack.x.push_back(pos_init.X());
    atrack.y.push_back(pos_init.Y());
    atrack.z.push_back(pos_init.Z());

    atrack.vx.push_back(dir_init.X());
    atrack.vy.push_back(dir_init.Y());
    atrack.vz.push_back(dir_init.Z());

    atrack.p.push_back(cal_loss.get_mom(T_init));
    // get charge ...
    double charge = 0;
    
    double dEdx = cal_loss.get_dEdx(T_init, step_size); // MeV
    T_init -= dEdx * step_size;
    charge = dEdx * step_size / (23.6 * units::eV) * lar.recombine_Birks(dEdx,88*units::kelvin,0.273*units::kilovolt/units::cm);
    // std::cout << T_init/units::MeV << " " << dEdx / (units::MeV/units::cm) << " " << lar.recombine_Birks(dEdx,88*units::kelvin,0.273*units::kilovolt/units::cm) << " " << charge << std::endl;
    atrack.Q.push_back(charge);
    pos_init.SetXYZ(pos_init.X() - dir_init.X() * step_size,
		    pos_init.Y() - dir_init.Y() * step_size,
		    pos_init.Z() - dir_init.Z() * step_size);
    atrack.N ++;
  }

  std::cout << atrack.N << std::endl;
}

void WireCell::MCS::MCSTrackSim(MCStrack& atrack, int particle_type, double T_init, TVector3 pos_init, TVector3 dir_init, double step_size){
  Eloss cal_loss(particle_type); // only for other particles ...

  LAr lar;

  // proper normalize ... 
  atrack.N = 0;
  dir_init *= 1./dir_init.Mag();

   while (T_init >0){
    atrack.x.push_back(pos_init.X());
    atrack.y.push_back(pos_init.Y());
    atrack.z.push_back(pos_init.Z());

    atrack.vx.push_back(dir_init.X());
    atrack.vy.push_back(dir_init.Y());
    atrack.vz.push_back(dir_init.Z());

    atrack.p.push_back(cal_loss.get_mom(T_init));
    // get charge ...
    double charge = 0;
    
    double dEdx = cal_loss.get_dEdx(T_init, step_size); // MeV

    
    
    charge = dEdx * step_size / (23.6 * units::eV) * lar.recombine_Birks(dEdx,88*units::kelvin,0.273*units::kilovolt/units::cm);

    // std::cout << "haha: " << dEdx / (units::MeV/units::cm) << " " << charge << std::endl;
    // std::cout << lar.recombine_Birks(1.8*units::MeV/units::cm,88*units::kelvin,0.273*units::kilovolt/units::cm);
    // std::cout << T_init/units::MeV << " " << dEdx / (units::MeV/units::cm) << " " << lar.recombine_Birks(dEdx,88*units::kelvin,0.273*units::kilovolt/units::cm) << " " << charge << std::endl;
    atrack.Q.push_back(charge);

    // four random number for multiple scattering ... 
    double z1= gRandom->Gaus(0,1);
    double z2= gRandom->Gaus(0,1);
    double z3= gRandom->Gaus(0,1);
    double z4= gRandom->Gaus(0,1);

    // average MCS angle ...
    double theta = cal_loss.get_mcs_angle(T_init, step_size);

    // std::cout << T_init/units::MeV << " " << theta << std::endl;
    
    double dz = step_size;
    double dx = z1* step_size * theta/std::sqrt(12) + z2* step_size* theta*0.5;
    double dy = z3 *step_size * theta/std::sqrt(12) + z4* step_size* theta*0.5;

    double thetaX = z2 * theta;
    double thetaY = z4 * theta;

    TVector3 v1(dx,dy,dz); // position ...
    TVector3 v2(tan(thetaX), tan(thetaY), 1); // angle ...

    RotateUz(dir_init,v1);
    pos_init.SetXYZ(pos_init.X() + v1.X(),
		    pos_init.Y() + v1.Y(),
		    pos_init.Z() + v1.Z());
    
    RotateUz(dir_init,v2);
    dir_init = v2.Unit();
    
    T_init -= dEdx * step_size;
        
    atrack.N ++;
  }

  std::cout << atrack.N << std::endl;
  
}

   
void WireCell::MCS::MultiScattSim(MCStrack& atrack, int N, std::vector<double> initpos, std::vector<double> initdir){
  atrack.clear();
  double x=0;
  double y=0;
  double z=0; // moving in z'-axis
  double thetaX =0;
  double thetaY =0;

  double stepLen = 2*units::cm;
  double varTheta = 3*units::mrad;
  
  double varThetaPlane = varTheta / std::sqrt(2);
  TVector3 direction;
  direction.SetMagThetaPhi(1, initdir[0], initdir[1]);
  while(atrack.N < N){
    double z1= gRandom->Gaus(0,1);
    double z2= gRandom->Gaus(0,1);
    double z3= gRandom->Gaus(0,1);
    double z4= gRandom->Gaus(0,1);

    double dx = stepLen *thetaX + z1*stepLen*varThetaPlane/std::sqrt(12) + z2* stepLen* varThetaPlane*0.5;
    double dy = stepLen *thetaY + z3*stepLen*varThetaPlane/std::sqrt(12) + z4* stepLen* varThetaPlane*0.5;
    z += stepLen;
    x += dx;
    y += dy;
    thetaX += z2 * varThetaPlane;
    thetaY += z4 * varThetaPlane;
    atrack.N ++;
    TVector3 v1(x,y,z);
    RotateUz(direction, v1);
    atrack.x.push_back(v1.X() + initpos[0]);
    atrack.y.push_back(v1.Y() + initpos[1]);
    atrack.z.push_back(v1.Z() + initpos[2]);

  }
}
