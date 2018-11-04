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

    double get_mean_ioniz_dEdx(double T);
    double get_mean_total_dEdx(double T);
    double get_MPV_dEdx(double T, double dx=0.3*units::cm);
    double get_mean_dEdx(double T, double tcut=0);
    
    
    double Density(double temperature);
    
    
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
