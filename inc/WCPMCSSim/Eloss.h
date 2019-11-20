#ifndef WIRECELL_ELOSS_H
#define WIRECELL_ELOSS_H

#include "WCPData/Units.h"

#include "TString.h"
#include "TGraph.h"
#include "TRandom.h"

namespace WCP {
  class Eloss {
  public:
    Eloss(int flag, TString filename = "input_data_files/proton_argon.dat");
    
    ~Eloss();

    // calculate mean values ... 
    double get_mean_ioniz_dEdx(double T);
    double get_mean_total_dEdx(double T);
    double get_MPV_dEdx(double T, double dx=0.3*units::cm);
    double get_mean_dEdx(double T, double tcut=0);

    // add random ...
    double get_dEdx(double T, double dx=0.3*units::cm);

    // calculate the multiple scattering angle ...
    double get_mcs_angle(double T, double dx=0.3*units::cm);
    
    // some unitility functions
    double Density(double temperature);
    double get_mom(double T);
    double get_kepa(double T, double dx = 0.3*units::cm);
    
  protected:
    int flag;
    double mass_p;
    double rho_lar;
    double fTemperature;
    double fZ, fA, fI, fSa, fSk, fSx0, fSx1, fScbar;
    double K, me;
    double fRadiationLength;
    
    double mass;
    double lifetime;

    std::vector<double> TE;
    std::vector<double> dEdx_rho;
    TGraph *g1; // ionization energy loss
    TGraph *g2; // total energy loss
  };
  
};

#endif
