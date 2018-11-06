#include "WireCellMCSSim/Eloss.h"

#include "TF1.h"

#include <fstream>
#include <iostream>

using namespace WireCell;

double WireCell::Eloss::Density(double temperature){
   // temperature is assumed to be in degrees Kelvin
    // density is nearly a linear function of temperature.  see the NIST tables for details
    // slope is between -6.2 and -6.1, intercept is 1928 kg/m^3
    // this parameterization will be good to better than 0.5%.
    // density is returned in g/cm^3
  double density = (-0.00615*temperature/units::kelvin + 1.928)*units::g/pow(units::cm,3);
  return density;
}



WireCell::Eloss::Eloss(int flag, TString filename)
  : flag(flag)
  , mass_p(938.272*units::MeV)
  , fTemperature(87.8*units::kelvin)
  , fZ(18)
  , fA(39.948*units::g/units::mole)
  , fI(188.0*units::eV)
  , fSa(0.1956)
  , fSk(3.0)
  , fSx0(0.2)
  , fSx1(3.0)
  , fScbar(5.2146)
  , K(0.307075*units::MeV*pow(units::cm,2)/units::mole)
  , me(0.510998918*units::MeV)
  , fRadiationLength(19.55*units::g/pow(units::cm,2)) 
{
  gRandom->SetSeed(0);
  
  rho_lar = Density(fTemperature);//1.396 * units::g/pow(units::cm,3);

  if (flag!=5){
    if (flag==1){ // pi+-
      mass = 139.570*units::MeV;
      lifetime = 2.6e-8*units::second; //sec
    }else if (flag==2){ // kaon +- 
      mass = 493.667*units::MeV;
      lifetime = 1.238e-8*units::second; //sec
    }else if (flag==3){ // muon+-
      mass = 105.658*units::MeV;
      lifetime = 2.2-6*units::second; // sec
    }else if (flag==4){
      mass = 938.272 * units::MeV; // proton
      lifetime = 100000.*units::second; //sec
    }
    
    
    
    std::ifstream infile(filename);
    double temp,beta;
    double temp_TE, temp_dEdx_rho;
    //while(!infile.eof()){
    for (Int_t i=0;i!=132;i++){
      infile >> temp_TE >> temp >> temp >> temp_dEdx_rho >> temp >> temp >> temp;
      temp_TE *= units::MeV;
      temp_dEdx_rho *= units::MeV/units::g * pow(units::cm,2);
      beta = sqrt(1-pow(mass_p/(mass_p+temp_TE),2));
      TE.push_back(mass/sqrt(1-beta*beta)-mass);
      dEdx_rho.push_back(temp_dEdx_rho*rho_lar);
    }
    Double_t x[10000],y[10000];
    for (size_t i=0;i+1 < TE.size();i++){
      x[i] = TE.at(i);
      y[i] = dEdx_rho.at(i);
    }
    int ncount = TE.size();
    //std::cout << ncount << std::endl;
    g1 = new TGraph(ncount,x,y);
    g2 = new TGraph(ncount,x,y);
  }else{
    //std::cout << " haha" << std::endl; 
    mass = me;
    lifetime = 100000.*units::second; //sec
    
    std::ifstream infile(filename);
    Double_t alpha_x[79],alpha_y[79],alpha_y1[79],temp;
    for (Int_t i=0;i!=79;i++){
      infile >> alpha_x[i] >> alpha_y[i] >> temp >> alpha_y1[i] >> temp;
      alpha_y[i] *= rho_lar * units::MeV/units::g * pow(units::cm,2);
      alpha_y1[i] *= rho_lar * units::MeV/units::g * pow(units::cm,2);
      alpha_x[i] *= units::MeV;
      TE.push_back(alpha_x[i]);
    }
    infile.close();
    
    g1 = new TGraph(79,alpha_x,alpha_y);
    g2 = new TGraph(79,alpha_x,alpha_y1);
  }
  
}

WireCell::Eloss::~Eloss(){
  delete g1;
  delete g2;
}

double WireCell::Eloss::get_mom(double T){
  double totE = T + mass;
  double mom = sqrt(pow(totE,2)-pow(mass,2));
  return mom;
}

double WireCell::Eloss::get_kepa(double T, double dx){
  double totE = T + mass;
  double mom = sqrt(pow(totE,2)-pow(mass,2));
  
  double bg = mom / mass;           // beta*gamma.
  double gamma = sqrt(1. + bg*bg);  // gamma.
  double beta = bg / gamma;         // beta (velocity).
  double mer = me / mass;   // electron mass / mass of incident particle.
  double tmax = 2.*me* bg*bg / (1. + 2.*gamma*mer + mer*mer);  // Maximum delta ray energy (MeV).
  
  double cosi = K/2. * fZ / fA /pow(beta,2) * dx * rho_lar;

  double kepa = cosi/tmax;
  return kepa;
}

double WireCell::Eloss::get_mcs_angle(double T, double dx){
  double totE = T + mass;
  double mom = sqrt(pow(totE,2)-pow(mass,2));
  
  double bg = mom / mass;           // beta*gamma.
  double gamma = sqrt(1. + bg*bg);  // gamma.
  double beta = bg / gamma;         // beta (velocity).

  double theta;
  theta = 13.6*units::MeV/beta/mom*sqrt(dx*rho_lar/fRadiationLength)*(1+0.038*log(dx*rho_lar/fRadiationLength/pow(beta,2)));
  return theta;
}


double WireCell::Eloss::get_dEdx(double T, double dx){
  double totE = T + mass;
  double mom = sqrt(pow(totE,2)-pow(mass,2));
  
  double bg = mom / mass;           // beta*gamma.
  double gamma = sqrt(1. + bg*bg);  // gamma.
  double beta = bg / gamma;         // beta (velocity).
  double mer = me / mass;   // electron mass / mass of incident particle.
  double tmax = 2.*me* bg*bg / (1. + 2.*gamma*mer + mer*mer);  // Maximum delta ray energy (MeV).
  
  double cosi = K/2. * fZ / fA /pow(beta,2) * dx * rho_lar;

  double kepa = cosi/tmax;


  
  //std::cout << cosi/units::MeV << " " << tmax/units::MeV << " " << kepa << std::endl;
  

  double dEdx=0;
  double gammap = 0.422784;
  //double gamma = 0.577215;
  double dEdx_ave = get_mean_dEdx(T); 

  if (kepa < 0.01){
    // Landau distribution
    double lambda = gRandom->Landau(0,1);
    dEdx = cosi*(lambda + log(kepa) + pow(beta,2)+gammap)/dx + dEdx_ave;
  }else if (kepa>=0.01 && kepa < 10){
    // Vavilov distribution
    TF1 f1("f1","TMath::Vavilov(x,[0],[1])",-5,40);
    f1.SetParameter(0,kepa);
    f1.SetParameter(1,beta*beta);
    double lambda = f1.GetRandom();
    dEdx = cosi*(lambda+log(kepa) + gammap + beta*beta)/dx + dEdx_ave;
  }else if (kepa >=10){
    // Gaussian distribution
    double sigma2 = cosi*tmax*(1-beta*beta/2.)/dx/dx;
    dEdx = gRandom->Gaus(dEdx_ave, sqrt(sigma2));
  }

  
  return dEdx;
}


double WireCell::Eloss::get_MPV_dEdx(double T, double dx){
  double totE = T + mass;
  double mom = sqrt(pow(totE,2)-pow(mass,2));
  
  double bg = mom / mass;           // beta*gamma.
  double gamma = sqrt(1. + bg*bg);  // gamma.
  double beta = bg / gamma;         // beta (velocity).

  // std::cout << gamma << " " << beta << " " << bg << std::endl;
  
  double psi = K/2.*fZ/fA/beta/beta*rho_lar*dx;
  
  double coef = log(2*me*bg*bg/fI) + log(psi/fI) + 0.2 - beta*beta;
  double x = log10(bg);
  double delta = 0.;
  if(x >= fSx0) {
    delta = 2. * log(10.) * x - fScbar;
    if(x < fSx1) {
      delta += fSa * pow(fSx1 - x, fSk);
    }
  }
  coef -= delta;
  
  return (psi * coef/dx);
}

double WireCell::Eloss::get_mean_dEdx(double T, double tcut){
  double totE = T + mass;
  double mom = sqrt(pow(totE,2)-pow(mass,2));

  // Calculate kinematic quantities.
  double bg = mom / mass;           // beta*gamma.
  double gamma = sqrt(1. + bg*bg);  // gamma.
  double beta = bg / gamma;         // beta (velocity).
  double mer = me / mass;   // electron mass / mass of incident particle.
  double tmax = 2.*me* bg*bg / (1. + 2.*gamma*mer + mer*mer);  // Maximum delta ray energy (MeV).
  // Make sure tcut does not exceed tmax.
  if(tcut == 0. || tcut > tmax) { tcut = tmax; }

  // Calculate density effect correction (delta).
  double x = log10(bg);
  double delta = 0.;
  if(x >= fSx0) {
    delta = 2. * log(10.) * x - fScbar;
    if(x < fSx1) {
      delta += fSa * pow(fSx1 - x, fSk);
    }
  }

  // Calculate stopping number.
  double B = 0.5 * log(2.*me*bg*bg*tcut / ( fI*fI))
    - 0.5 * beta*beta * (1. + tcut / tmax) - 0.5 * delta;
  // Don't let the stopping number become negative.
  if(B < 1.) B = 1.;
  
  // Calculate dE/dx.
  double dedx = rho_lar * K*fZ*B / (fA * beta*beta);
  return dedx;   // MeV/cm.
}


double WireCell::Eloss::get_mean_ioniz_dEdx(double T){
  if (T > TE.at(TE.size()-1))
    T = TE.at(TE.size()-1);
  return g1->Eval(T);
}

double WireCell::Eloss::get_mean_total_dEdx(double T){
  if (T > TE.at(TE.size()-1))
    T = TE.at(TE.size()-1);
  return g2->Eval(T);
}


