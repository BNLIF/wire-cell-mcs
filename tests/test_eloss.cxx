#include "WireCellMCSSim/Eloss.h"
#include "WireCellData/LAr.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"


#include <iostream>

using namespace WireCell;
using namespace std;




int main(int argc, char* argv[])
{
  Eloss cal_eloss(3); // pion, kaon, muon, proton, electron

  std::cout << cal_eloss.get_mean_ioniz_dEdx(100*units::MeV) / units::MeV * units::cm << " MeV/cm " << cal_eloss.get_MPV_dEdx(100*units::MeV, 0.03*units::cm)/ units::MeV * units::cm << " " << cal_eloss.get_mean_dEdx(100*units::MeV, 0)/units::MeV * units::cm << std::endl;

  std::cout << cal_eloss.get_dEdx(10*units::MeV,0.03*units::cm)/ units::MeV * units::cm << std::endl;
  
  std::cout << "MCS angle: " << cal_eloss.get_mcs_angle(10*units::MeV,0.1*units::cm) << std::endl;
  
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

  // simulate eloss ...
  TFile *file = new TFile("eloss.root","RECREATE");
  TTree *T_ele = new TTree("T_ele","T_ele");

  TTree *T_pi = new TTree("T_pi","T_pi");
  TTree *T_mu = new TTree("T_mu","T_mu");
  TTree *T_p = new TTree("T_p","T_p");
  TTree *T_k = new TTree("T_k","T_k");

  TH1F *h1 = new TH1F("h1","h1",2000,0,20);
  h1->SetDirectory(file);
  for (int i=0;i!=10000;i++){
    h1->Fill( cal_eloss.get_dEdx(100*units::MeV,30*units::cm)/ units::MeV * units::cm);
  }
  
  Double_t Ekin, Lstop, dEodx;
  Double_t dEodx_MPV;
  T_ele->Branch("T",&Ekin,"T/D");
  T_ele->Branch("L",&Lstop,"L/D");
  T_ele->Branch("dEodx",&dEodx,"dEodx/D");
  T_ele->Branch("dEodx_MPV",&dEodx_MPV,"dEodx_MPV/D");
  T_ele->SetDirectory(file);

  T_pi->Branch("T",&Ekin,"T/D");
  T_pi->Branch("L",&Lstop,"L/D");
  T_pi->Branch("dEodx",&dEodx,"dEodx/D");
  T_pi->Branch("dEodx_MPV",&dEodx_MPV,"dEodx_MPV/D");
  T_pi->SetDirectory(file);
  
  T_p->Branch("T",&Ekin,"T/D");
  T_p->Branch("L",&Lstop,"L/D");
  T_p->Branch("dEodx",&dEodx,"dEodx/D");
  T_p->Branch("dEodx_MPV",&dEodx_MPV,"dEodx_MPV/D");
  T_p->SetDirectory(file);

  T_k->Branch("T",&Ekin,"T/D");
  T_k->Branch("L",&Lstop,"L/D");
  T_k->Branch("dEodx",&dEodx,"dEodx/D");
  T_k->Branch("dEodx_MPV",&dEodx_MPV,"dEodx_MPV/D");
  T_k->SetDirectory(file);

  T_mu->Branch("T",&Ekin,"T/D");
  T_mu->Branch("L",&Lstop,"L/D");
  T_mu->Branch("dEodx",&dEodx,"dEodx/D");
  T_mu->Branch("dEodx_MPV",&dEodx_MPV,"dEodx_MPV/D");
  T_mu->SetDirectory(file);
  
  Eloss ele_eloss(5,"input_data_files/electron_argon.dat");

  std::cout << ele_eloss.get_mean_ioniz_dEdx(60*units::MeV) / units::MeV * units::cm << " MeV/cm" << std::endl;
  
  // electron simulation ...
  Ekin = 60; // MeV
  Lstop = 0;
  for(int i=0;i!=60*50000;i++){
    Double_t dE = 1./50000.; // MeV
    
    dEodx = ele_eloss.get_mean_dEdx(Ekin*units::MeV) / units::MeV * units::cm ;
    dEodx_MPV = ele_eloss.get_MPV_dEdx(Ekin*units::MeV,3*units::mm)/ units::MeV * units::cm ;
    
    double dEodx_all = ele_eloss.get_mean_total_dEdx(Ekin*units::MeV) / units::MeV * units::cm;
    Lstop += dE/dEodx_all ;
    if (i%5000==0)
    //    if ((int(Ekin/0.01))%10==0)
      T_ele->Fill();
    Ekin -= dE; // MeV ...
  }

  Eloss pion_eloss(1); // pion, kaon, muon, proton, electron
  Eloss kaon_eloss(2); // pion, kaon, muon, proton, electron
  Eloss muon_eloss(3); // pion, kaon, muon, proton, electron
  Eloss prot_eloss(4); // pion, kaon, muon, proton, electron

  Ekin = 200;
    Lstop = 0;
  for(int i=0;i!=200*5000;i++){
    Double_t dE = 1./5000.; // MeV
    dEodx = prot_eloss.get_mean_dEdx(Ekin*units::MeV) / units::MeV * units::cm ;
    dEodx_MPV = prot_eloss.get_MPV_dEdx(Ekin*units::MeV,3*units::mm)/ units::MeV * units::cm ;
    double dEodx_all = prot_eloss.get_mean_total_dEdx(Ekin*units::MeV) / units::MeV * units::cm;
    Lstop += dE/dEodx_all ;
    if (i%500==0)
    //    if ((int(Ekin/0.01))%10==0)
      T_p->Fill();
    Ekin -= dE; // MeV ...
  }

  Ekin = 200;
    Lstop = 0;
  for(int i=0;i!=200*5000;i++){
    Double_t dE = 1./5000.; // MeV
    dEodx = kaon_eloss.get_mean_dEdx(Ekin*units::MeV) / units::MeV * units::cm ;
    dEodx_MPV = kaon_eloss.get_MPV_dEdx(Ekin*units::MeV,3*units::mm)/ units::MeV * units::cm ;
    double dEodx_all = kaon_eloss.get_mean_total_dEdx(Ekin*units::MeV) / units::MeV * units::cm;
    Lstop += dE/dEodx_all ;
    if (i%500==0)
    //    if ((int(Ekin/0.01))%10==0)
      T_k->Fill();
    Ekin -= dE; // MeV ...
  }

  Ekin = 200;
    Lstop = 0;
  for(int i=0;i!=200*5000;i++){
    Double_t dE = 1./5000.; // MeV
    dEodx = pion_eloss.get_mean_dEdx(Ekin*units::MeV) / units::MeV * units::cm ;
    dEodx_MPV = pion_eloss.get_MPV_dEdx(Ekin*units::MeV,3*units::mm)/ units::MeV * units::cm ;
    double dEodx_all = pion_eloss.get_mean_total_dEdx(Ekin*units::MeV) / units::MeV * units::cm;
    Lstop += dE/dEodx_all ;
    if (i%500==0)
    //    if ((int(Ekin/0.01))%10==0)
      T_pi->Fill();
    Ekin -= dE; // MeV ...
  }

  Ekin = 200;
    Lstop = 0;
  for(int i=0;i!=200*5000;i++){
    Double_t dE = 1./5000.; // MeV
    dEodx = muon_eloss.get_mean_dEdx(Ekin*units::MeV) / units::MeV * units::cm ;
    dEodx_MPV = muon_eloss.get_MPV_dEdx(Ekin*units::MeV,3*units::mm)/ units::MeV * units::cm ;
    double dEodx_all = muon_eloss.get_mean_total_dEdx(Ekin*units::MeV) / units::MeV * units::cm;
    Lstop += dE/dEodx_all ;
    if (i%500==0)
    //    if ((int(Ekin/0.01))%10==0)
      T_mu->Fill();
    Ekin -= dE; // MeV ...
  }
  
  
  file->Write();
  file->Close();
  
  return 0;
  //  std::cout << "test" << std::endl;
}
