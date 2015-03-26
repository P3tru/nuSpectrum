#include <stdio.h>
#include <iostream>

#include "TFile.h"
#include "TH1D.h"

#include "utils.hh"
#include "classNuSpectrum.hh"

// ************** Compute the Neutrinos flux in the Beam from ILL reactor ************** //
double neutrino_flux(double E_ve)
{
  double E_235U=201.9; // 201.9 + 0.5 MeV arxiv:hep-ph/0410100
  double Pth=3.557658626964e20; // Power Thermal Spectrum 57MW at ILL converted in eV
  double S_235U=1; // Reference antineutrino spectrum
  double alpha[6]; alpha[0]=3.217;alpha[1]=-3.111;alpha[2]=1.395;alpha[3]=-0.3690;alpha[4]=0.04445;alpha[5]=-0.002053; // arxiv:hep-ex/1101.2663
  
  for(int i=0;i<6;i++){ S_235U*=exp(alpha[i]*pow(E_ve,i));}
  
  double flux=Pth*S_235U/E_235U;
  //  double flux=S_235U;
  //  printf("Nu flux = %e \n", flux);
  return flux;
    
}

// ************** Compute antineutrino cross-section  ************** //
double neutrino_sigma(double E_ve)
{
  double K=0.961e-43; // A.Pichlmaier et al Phys Lett B693 221 (2010)
  double Delta=1.293; // Mn=939,565 MeV Mp=938,272 MeV
  double Me = 0.511; // 511 KeV

  double cross = K*(E_ve-Delta)*sqrt((E_ve-Delta)*(E_ve-Delta)-Me*Me);
  //  double cross = 6.69e-43; // Mean value for cut at 2 MeV
  return cross;
}


int main(int argc, char **argv){

  NuSpectrum spectrum;

  int Emin = 2; //MeV
  int Emax = 8; //MeV
  int ERes = 10; //%
  int nbBins=(Emax-Emin)*ERes;
  double intBins=1/(double)ERes;
  double E[nbBins];
  for(int i=0;i<nbBins;i++) E[i]=Emin+intBins*i;
  
  TH1D *hCell1 = new TH1D("hCell1","hCell1",nbBins+1,Emin,Emax);
  TH1D *hCell2 = new TH1D("hCell2","hCell2",nbBins+1,Emin,Emax);
  TH1D *hCell3 = new TH1D("hCell3","hCell3",nbBins+1,Emin,Emax);
  TH1D *hCell4 = new TH1D("hCell4","hCell4",nbBins+1,Emin,Emax);
  TH1D *hCell5 = new TH1D("hCell5","hCell5",nbBins+1,Emin,Emax);
  TH1D *hCell6 = new TH1D("hCell6","hCell6",nbBins+1,Emin,Emax);
  
  double posCell1=9; // position Cell 1 9m
  double integral1;
  double posCell2=posCell1+0.4;
  double integral2;
  double posCell3=posCell1+0.8;
  double integral3;
  double posCell4=posCell1+1.2;
  double integral4;
  double posCell5=posCell1+1.6;
  double integral5;
  double posCell6=posCell1+2.0;
  double integral6;
  
  const int n_int=1000;
  const double Np=6.3e28;
  double z[n_int]; double w[n_int];
  
  for(int i=0;i<nbBins;i++){
    gauleg(E[i],E[i]+intBins,n_int,z,w);	  
    integral1=0;
    integral2=0;
    integral3=0;
    integral4=0;
    integral5=0;
    integral6=0;
    for (int j=0;j<n_int;j++){
      integral1+=w[j]*Np*neutrino_flux(z[j])*neutrino_sigma(z[j])*spectrum.Osc4Nu(z[j],posCell1);
      integral2+=w[j]*Np*neutrino_flux(z[j])*neutrino_sigma(z[j])*spectrum.Osc4Nu(z[j],posCell2);
      integral3+=w[j]*Np*neutrino_flux(z[j])*neutrino_sigma(z[j])*spectrum.Osc4Nu(z[j],posCell3);
      integral4+=w[j]*Np*neutrino_flux(z[j])*neutrino_sigma(z[j])*spectrum.Osc4Nu(z[j],posCell4);
      integral5+=w[j]*Np*neutrino_flux(z[j])*neutrino_sigma(z[j])*spectrum.Osc4Nu(z[j],posCell5);
      integral6+=w[j]*Np*neutrino_flux(z[j])*neutrino_sigma(z[j])*spectrum.Osc4Nu(z[j],posCell6);
    }
    hCell1->Fill(E[i],integral1/intBins);
    hCell2->Fill(E[i],integral2/intBins);
    hCell3->Fill(E[i],integral3/intBins);
    hCell4->Fill(E[i],integral4/intBins);
    hCell5->Fill(E[i],integral5/intBins);
    hCell6->Fill(E[i],integral6/intBins);
    printf("E=%f int1=%e int2=%f int3=%f int4=%f int5=%f int6=%f\n",
     	   E[i],integral1/intBins,integral2/intBins,integral3/intBins,integral4/intBins,integral5/intBins,integral6/intBins);
  }

  for(int i=0;i<nbBins*2;i++){
    hCell1->SetBinError(i,0);
    hCell2->SetBinError(i,0);
    hCell3->SetBinError(i,0);
    hCell4->SetBinError(i,0);
    hCell5->SetBinError(i,0);
    hCell6->SetBinError(i,0);
  }

  TFile f("simuCells.root","RECREATE");
  hCell1->Write();
  hCell2->Write();
  hCell3->Write();
  hCell4->Write();
  hCell5->Write();
  hCell6->Write();
  
  return 0;
}
