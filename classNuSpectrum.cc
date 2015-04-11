#include "classNuSpectrum.hh"
#include "utils.hh"

namespace std {} using namespace std;

NuSpectrum::NuSpectrum()
{
  // Mixing Angles
  s12=sqrt(0.306); c12=sqrt(1-0.306);	
  s13=sin(asin(sqrt(0.092))/2); c13=sqrt(1-s13*s13);	
  s23=sqrt(0.42); c23=sqrt(1-0.42);	
  s34=0; c34=1;
  s24=0; c24=1;
  s14=sin(asin(sqrt(0.17))/2); c14=sqrt(1-s14*s14);

  //Mixing matrix
  mat13.ResizeTo(4,4);
  mat12.ResizeTo(4,4);
  mat23.ResizeTo(4,4);
  mat34.ResizeTo(4,4);
  mat24.ResizeTo(4,4);
  mat14.ResizeTo(4,4);
  matU.ResizeTo(4,4);
  //Fill mixing matrix
  mat13(0,0)=c13; mat13(0,1)=0; mat13(0,2)=s13; mat13(0,3)=0;
  mat13(1,0)=0; mat13(1,1)=1; mat13(1,2)=0; mat13(1,3)=0;
  mat13(2,0)=-s13; mat13(2,1)=0; mat13(2,2)=c13; mat13(2,3)=0;
  mat13(3,0)=0; mat13(3,1)=0; mat13(3,2)=0; mat13(3,3)=1;

  mat12(0,0)=c12; mat12(0,1)=s12; mat12(0,2)=0; mat12(0,3)=0;
  mat12(1,0)=-s12; mat12(1,1)=c12; mat12(1,2)=0; mat12(1,3)=0;		
  mat12(2,0)=0; mat12(2,1)=0; mat12(2,2)=1; mat12(2,3)=0;		
  mat12(3,0)=0; mat12(3,1)=0; mat12(3,2)=0; mat12(3,3)=1;	

  mat23(0,0)=1; mat23(0,1)=0; mat23(0,2)=0; mat23(0,3)=0;
  mat23(1,0)=0; mat23(1,1)=c23; mat23(1,2)=s23; mat23(1,3)=0;
  mat23(2,0)=0; mat23(2,1)=-s23; mat23(2,2)=c23; mat23(2,3)=0;
  mat23(3,0)=0; mat23(3,1)=0; mat23(3,2)=0; mat23(3,3)=1;	

  mat34(0,0)=1; mat34(0,1)=0; mat34(0,2)=0; mat34(0,3)=0;
  mat34(1,0)=0; mat34(1,1)=1; mat34(1,2)=0; mat34(1,3)=0;
  mat34(2,0)=0; mat34(2,1)=0; mat34(2,2)=c34; mat34(2,3)=s34;
  mat34(3,0)=0; mat34(3,1)=0; mat34(3,2)=-s34; mat34(3,3)=c34;

  mat24(0,0)=1; mat24(0,1)=0;mat24(0,2)=0;mat24(0,3)=0;
  mat24(1,0)=0; mat24(1,1)=c24; mat24(1,2)=0; mat24(1,3)=s24;
  mat24(2,0)=0; mat24(2,1)=0;mat24(2,2)=1;mat24(2,3)=0;
  mat24(3,0)=0; mat24(3,1)=-s24; mat24(3,2)=0; mat24(3,3)=c24;

  mat14(0,0)=c14; mat14(0,1)=0; mat14(0,2)=0; mat14(0,3)=s14;
  mat14(1,0)=0; mat14(1,1)=1; mat14(1,2)=0; mat14(1,3)=0;
  mat14(2,0)=0; mat14(2,1)=0; mat14(2,2)=1; mat14(2,3)=0; 
  mat14(3,0)=-s14; mat14(3,1)=0; mat14(3,2)=0; mat14(3,3)=c14;

  mixMatrix.ResizeTo(4,4);
  mixMatrix=mat34*mat24*mat14*mat23*mat13*mat12;

  // Delta Masses
  deltam21=7.58e-5;
  deltam31=2.35e-5;
  deltam32=2.35e-3;
  deltam41=2.3;
  deltam42=0;
  deltam43=0;

  for(int i=0;i<4;i++) for(int j=0;j<4;j++) deltam[i][j]=0;
  deltam[1][0]=deltam21;
  deltam[2][1]=deltam32;
  deltam[2][0]=deltam31;
  deltam[3][2]=deltam43;
  deltam[3][1]=deltam42;
  deltam[3][0]=deltam41;
}

NuSpectrum::~NuSpectrum(){
}

double NuSpectrum::Osc3Nu(double EE, double LL){
  double sin2theta[3];
  for(int i=0;i<3;i++) sin2theta[i]=0;

  for (int i=1;i<3;i++){
    for (int j=0;j<i;j++){
      sin2theta[i]+=4*mixMatrix(0,i)*mixMatrix(0,i)*mixMatrix(0,j)*mixMatrix(0,j)*pow(sin(1.27*deltam[i][j]*LL/EE),2);
    }
  }  

  return 1-sin2theta[1]-sin2theta[2];
}


double NuSpectrum::Osc4Nu(double EE, double LL){
  double sin2theta[4];
  for(int i=0;i<4;i++) sin2theta[i]=0;
	for (int i=1;i<4;i++){
		for (int j=0;j<i;j++){
		sin2theta[i]+=4*mixMatrix(0,i)*mixMatrix(0,i)*mixMatrix(0,j)*mixMatrix(0,j)*pow(sin(1.27*deltam[i][j]*LL/EE),2);
		}
	}

	double f=1-sin2theta[1]-sin2theta[2]-sin2theta[3];
	return f;
}

TGraph NuSpectrum::IntGraph(double minE, double maxE, int nbNu){
  if(nbNu == -1) return -1;
  // const long int n=100000; // Domain of the Neutrino Spectrum in Meter
  const int n=10000; // Domain of the Neutrino Spectrum in Meter
  const int n_int=1000; // Number of integrations point
  double x[2*n],y[2*n];
  for(int i=0;i<n;i++){
    x[i]=i;
    y[i]=0;
    // printf("x[%d]=%.0f\n",i,x[i]);
  }
  for(int i=0;i<n;i++){
    x[n+i]=n*(0.01*i+1);
    y[n+i]=0;
    // printf("x[%d]=%.0f\n",n+i,x[n+i]);
  }

  double integral;
  double z[n_int]; double w[n_int];
  gauleg(minE,maxE,n_int,z,w);
  for (int i=0;i<2*n;i++){
    integral=0;
    for (int j=0;j<n_int;j++){
      if(nbNu == 3) integral+=w[j]*Osc3Nu(z[j],x[i]);
      if(nbNu == 4) integral+=w[j]*Osc4Nu(z[j],x[i]);
    }
    y[i]=integral/(maxE-minE);
  }  
  TGraph gr(2*n,x,y);
  gr.SetLineColor(2);
  gr.SetLineWidth(2);
  gr.SetMarkerColor(2);
  gr.SetMarkerStyle(7);
  gr.SetTitle("Neutrino Oscillation Spectrum");
  gr.GetXaxis()->SetTitle("Distance from reactor (m)");
  gr.GetYaxis()->SetTitle("Surviving Probability Ve -> Ve");
  gr.GetYaxis()->SetRange(0,1);
  return gr;
}

TGraph NuSpectrum::Graph(double E, int nbNu){
  if(nbNu == -1) return -1;
  const int n=100000; // Domain of the Neutrino Spectrum in Meter
  double x[n],y[n];
  for(int i=0;i<n;i++){
    x[i]=i;
    if(nbNu == 3) y[i]=Osc3Nu(E,x[i]);
    if(nbNu == 4) y[i]=Osc4Nu(E,x[i]);
  }
  TGraph gr(n,x,y);
  return gr;
}

void NuSpectrum::SetOptionsGr(TGraph *gr){
  gr->SetTitle("Neutrino Oscillation Spectrum");
  gr->GetXaxis()->SetTitle("Distance from reactor (m)");
  gr->GetYaxis()->SetTitle("Surviving Probability Ve -> Ve");
  gr->GetYaxis()->SetRange(0,1);
}

void NuSpectrum::simuCells(NuSpectrum spectrum, TFile *f)
{
	int Emin = 2; //MeV
	int Emax = 8; //MeV
	int ERes = 10; //%
	int nbBins=(Emin-Emax)*ERes;
	double intBins=1/(double)ERes;
	double E[nbBins];
	for(int i=0;i<nbBins;i++) E[i]=Emin+intBins*i;

	TH1D *hCell1 = new TH1D("hCell1","hCell1",nbBins,Emin,Emax);
	TH1D *hCell2 = new TH1D("hCell2","hCell2",nbBins,Emin,Emax);
	TH1D *hCell3 = new TH1D("hCell3","hCell3",nbBins,Emin,Emax);
	TH1D *hCell4 = new TH1D("hCell4","hCell4",nbBins,Emin,Emax);
	TH1D *hCell5 = new TH1D("hCell5","hCell5",nbBins,Emin,Emax);
	TH1D *hCell6 = new TH1D("hCell6","hCell6",nbBins,Emin,Emax);

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

	const int n_int=100;
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
	    integral1+=w[j]*spectrum.Osc4Nu(z[j],posCell1);
	    integral2+=w[j]*spectrum.Osc4Nu(z[j],posCell2);
	    integral3+=w[j]*spectrum.Osc4Nu(z[j],posCell3);
	    integral4+=w[j]*spectrum.Osc4Nu(z[j],posCell4);
	    integral5+=w[j]*spectrum.Osc4Nu(z[j],posCell5);
	    integral6+=w[j]*spectrum.Osc4Nu(z[j],posCell6);
	  }
	  hCell1->Fill(E[i],integral1/intBins);
	  hCell2->Fill(E[i],integral2/intBins);
	  hCell3->Fill(E[i],integral3/intBins);
	  hCell4->Fill(E[i],integral4/intBins);
	  hCell5->Fill(E[i],integral5/intBins);
	  hCell6->Fill(E[i],integral6/intBins);
	}

	f = new TFile("simuCells.root","RECREATE");
        hCell1->Write();
        hCell2->Write();
        hCell3->Write();
        hCell4->Write();
        hCell5->Write();
        hCell6->Write();
} // end simuCells


