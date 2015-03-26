#ifndef NuSpectrum_h
#define NuSpectrum_h

//User Includes
#include <stdio.h>
#include <iostream>
#include <math.h>

#include <TGraph.h>
#include <TMath.h> 
#include <TMatrixD.h>
#include <TAxis.h>
#include <TFile.h>
#include <TH1D.h>

class NuSpectrum{  
private:	//Data members
  
  // Mixing Angles
  double s12,c12;
  double s13,c13;
  double s23,c23;
  double s14,c14;
  double s24,c24;
  double s34,c34;
  // Mixing Matrices
  TMatrixD mat23;
  TMatrixD mat13;
  TMatrixD mat12;
  TMatrixD matU;
  TMatrixD mat34;
  TMatrixD mat24;
  TMatrixD mat14;
  //Define Mixing Matrix
  TMatrixD mixMatrix;

  // Delta Masses
  double deltam21;
  double deltam31;
  double deltam32;
  double deltam41;
  double deltam42;
  double deltam43;
  //Define Mixing Matrix
  double deltam[4][4];

public:		//Methods
  NuSpectrum();  //Constructor
  ~NuSpectrum(); //Destructor
  
  // Compute Ve->Ve surviving probability in 3-neutrinos model
  double Osc3Nu(double EE, double LL);

  // Compute Ve->Ve surviving probability in 4-neutrinos model
  double Osc4Nu(double EE, double LL);

  // Return integrated spectrum w.r.t distance between minE and maxE
  TGraph IntGraph(double minE=2, double maxE=6, int nbNu=-1);

  // Return integrated spectrum w.r.t distance between minE and maxE
  TGraph Graph(double E=1, int nbNu=-1);

  // Set Options for Graphs
  void SetOptionsGr(TGraph *gr);

  void simuCells(NuSpectrum spectrum, TFile *f);

};

#endif
