#include "WireCellMCSSim/Eloss.h"

#include <fstream>

using namespace WireCell;

WireCell::Eloss::Eloss(int flag, TString filename){
  mass_p = 938.272*units::MeV;
  rho_lar = 1.396 * units::g/pow(units::cm,3);

  if (flag==1){ // pi+-
    mass_pi = 139.570*units::MeV;
    pi_lifetime = 2.6e-8*units::second; //sec
  }else if (flag==2){ // kaon +- 
    mass_pi = 493.667*units::MeV;
    pi_lifetime = 1.238e-8*units::second; //sec
  }else if (flag==3){ // muon+-
    mass_pi = 105.658*units::MeV;
    pi_lifetime = 2.2-6*units::second; // sec
  }else if (flag==4){
    mass_pi = 938.272 * units::MeV; // proton
    pi_lifetime = 100000.*units::second; //sec
  }
  

  
  std::ifstream infile(filename);
  double temp,beta;
  double temp_TE, temp_dEdx_rho;
  while(!infile.eof()){
  //for (Int_t i=0;i!=132;i++){
    infile >> temp_TE >> temp >> temp >> temp_dEdx_rho >> temp >> temp >> temp;
    temp_TE *= units::MeV;
    //    temp_dEdx_rho = ;
    beta = sqrt(1-pow(mass_p/(mass_p+temp_TE),2));
    TE.push_back(mass_pi/sqrt(1-beta*beta)-mass_pi);
    dEdx_rho.push_back(temp_dEdx_rho*rho_lar);
  }
  Double_t x[10000],y[10000];
  for (size_t i=0;i+1 < TE.size();i++){
    x[i] = TE.at(i);
    y[i] = dEdx_rho.at(i);
  }
  int ncount = TE.size()-1;
  g1 = new TGraph(ncount,x,y);
  
}

WireCell::Eloss::~Eloss(){
  delete g1;
}

double WireCell::Eloss::get_mean_dEdx(double T){
  if (T > TE.at(TE.size()-1))
    T = TE.at(TE.size()-1);
  return g1->Eval(T);
}
