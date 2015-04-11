#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <signal.h>
#include <math.h>

#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TApplication.h>

#include "utils.hh"
#include "classNuSpectrum.hh"

void Interrupt(int arg){

  printf("got a control-C, exiting\n"); exit(0);
  
}

void usage(char *name){

  printf("\nUsage: %s -h[elp] file :\n",name);
  printf("\t-h: this message\n");
  printf("\n");
  
}

void simulationCells(NuSpectrum spectrum){
    int Emin = 2; //MeV
    int Emax = 8; //MeV
    int ERes = 10; //%
    int nbBins = 0;
    double intBins = 0;
    if(Emin && Emax) nbBins=(Emax-Emin)*ERes;
    if(ERes !=0 ) intBins=1/(double)ERes;
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
}

int main(int argc, char **argv){

  // signal(SIGINT,Interrupt);
  // argc=0;
  // TApplication theApp("App",&argc,argv);
  // theApp.Run();
  
  bool SimulationCells = false;  // Simulation des spectres en E dans les 6 cellules

  NuSpectrum spectrum;
  TGraph gr;
  gr = spectrum.IntGraph(2,8,4);

  TFile f("test.root","RECREATE");
  gr.Write();
  
  if(SimulationCells){ // Simulation des spectres en E dans les 6 cellules
    simulationCells(spectrum);
  }

  return 0;
}
