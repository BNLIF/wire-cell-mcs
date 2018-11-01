#ifndef WIRECELL_ELOSS_H
#define WIRECELL_ELOSS_H

#include "WireCellData/Units.h"

#include "TString.h"
#include "TGraph.h"

namespace WireCell {
  class Eloss {
  public:
    Eloss(int flag, TString filename = "input_data_files/proton_argon.dat");
    
    ~Eloss();

  protected:
    double mass_p;
    double rho_lar;
    
    double mass_pi;
    double pi_lifetime;

    std::vector<double> TE;
    std::vector<double> dEdx_rho;
    TGraph *g1;
  };
  
};

#endif
