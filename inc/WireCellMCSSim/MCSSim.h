#ifndef WIRECELL_MCSSIM_H
#define WIRECELL_MCSSIM_H

#include "WireCellData/Units.h"
#include "TVector3.h"
#include "TRandom.h"

#include <vector>

namespace WireCell {

  struct MCStrack{
    int N; // number of vertices 
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    
    std::vector<double> vx; // momentum direction
    std::vector<double> vy; // momentum direction
    std::vector<double> vz; // momentum direction
    std::vector<double> p; // momentum
    std::vector<double> Q; // charge deposition

    double mass; // mass of the particle
    double charge; // charge of the particle 
    
    void clear(){
      N = 0;
      x.clear();
      y.clear();
      z.clear();
    }
  };

  
  
  void RotateUz(TVector3& direction, TVector3& v1);

  void MultiScattSim(MCStrack& atrack, int N, std::vector<double> initpos, std::vector<double> initdir);

  void LineTrackSim(MCStrack& atrack, int particle_type, double T_init, TVector3 pos_init, TVector3 dir_init, double step_size);
  
};

#endif
