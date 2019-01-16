// Created on 10/31/2018
// by Wenqiang Gu (wgu@bnl.gov)
// to compile it standalone: g++ Mcs.cxx -o Mcs.exe `root-config --cflags --libs`

#include <iostream>
#include <vector>
#include <cmath>
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include "TH1F.h"

namespace MCS{
  double mrad = 1E-3;
  double cm = 1E-2;
  double m = 1.;
  double stepLen = 0.3*cm;
  double varTheta = 3*mrad;

  struct track{
    int N; // number of vertices 
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    void clear(){
      N = 0;
      x.clear(); y.clear(); z.clear();
    }
  };
};

void RotateUz(TVector3 direction, TVector3& v1){
  TVector3 unit = direction.Unit();
  v1.RotateUz(unit); 
}

// c.f. http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf
// vertex and direction change in the local coordinate (moving in z-axis)
void McsLocal(std::vector<double>& lMcs){
    double z1= gRandom->Gaus(0,1);
    double z2= gRandom->Gaus(0,1);
    double z3= gRandom->Gaus(0,1);
    double z4= gRandom->Gaus(0,1);
	double varThetaPlane = MCS::varTheta / std::sqrt(2);
    double dx = z1*MCS::stepLen*varThetaPlane/std::sqrt(12) + z2* MCS::stepLen* varThetaPlane*0.5;
    double dy = z3*MCS::stepLen*varThetaPlane/std::sqrt(12) + z4* MCS::stepLen* varThetaPlane*0.5;
    double dz = MCS::stepLen;
    double thx = z2 * varThetaPlane; // theta in x-z plane
    double thy = z4 * varThetaPlane;
    // double tan2th = pow(tan(thx),2) + pow(tan(thy),2);
    // double cos2th = 1./(1 + tan2th);
    // double costh = sqrt(cos2th);
    lMcs[0] = dx;
    lMcs[1] = dy;
    lMcs[2] = dz;
    lMcs[3] = tan(thx); //costh*tan(thx);
    lMcs[4] = tan(thy); //costh*tan(thy);
    lMcs[5] = 1.0; //costh;

    // return std::vector<double> v{dx,dy,dz, costh*tan(thx), costh*tan(thy), costh};// position and direction after mcs (in local corordinate)
}

void McsGlobal(MCS::track& atrack, int N, std::vector<double> initpos, std::vector<double> initdir, TH1F* h){
  atrack.clear();
  
  TVector3 prevXyz, prevDir;
  prevXyz.SetXYZ(initpos[0], initpos[1], initpos[2]);
  prevDir.SetMagThetaPhi(initdir[0], initdir[1], initdir[2]);
  while(atrack.N < N){

    std::vector<double> localMcs{0,0,0,0,0,0};
    McsLocal(localMcs);
    TVector3 dXyz(localMcs[0], localMcs[1], localMcs[2]);
    TVector3 Dir(localMcs[3], localMcs[4], localMcs[5]);

    RotateUz(prevDir, dXyz); // local dr (dXyz) to global dr
    RotateUz(prevDir, Dir); // local direction (Dir) to global direction

h->Fill(dXyz.Angle(prevDir));

    atrack.x.push_back(dXyz[0] + prevXyz[0]);
    atrack.y.push_back(dXyz[1] + prevXyz[1]);
    atrack.z.push_back(dXyz[2] + prevXyz[2]);

    atrack.N ++;
    prevDir = Dir;
    prevXyz += dXyz;
  }
}



int main(){
  MCS::track atrack;

  TFile* ofile = new TFile("mcs-tracks.root","RECREATE");
  TTree* T = new TTree("T","tracks and vertices");
  T->Branch("N", &atrack.N, "N/I");
  T->Branch("x", &atrack.x);
  T->Branch("y", &atrack.y);
  T->Branch("z", &atrack.z);

TH1F* h = new TH1F("h","",100,-0.01,0.01);
  for(int i=0; i<1000/*ntracks*/; i++){
    double x = gRandom->Uniform(2.56*MCS::m);
    double y = gRandom->Uniform(2.3*MCS::m);
    double z = gRandom->Uniform(10.4*MCS::m);
    double cosTheta = gRandom->Uniform(-1,1);
    double phi = gRandom->Uniform(0,2*3.1415926);
    std::vector<double> initpos{x,y,z};
    std::vector<double> initdir{1, acos(cosTheta), phi};

    McsGlobal(atrack, 100/*nvertices*/, initpos, initdir, h);
    T->Fill();
  }
  T->Write();
h->Write();
  ofile->Close();
}
