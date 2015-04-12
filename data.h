#ifndef DATA_H
#define DATA_H

#include "TGraphErrors.h"

TGraphErrors *pauloVerdeData(); // 820 metre. 2001
TGraphErrors *choozData(); // 1050 metre.
TGraphErrors *SRPD96ata(); // 18 et 24. metre.
TGraphErrors *krasnoyarskData(); // 33., 92. m, 1987 et 57. m 1994
TGraphErrors *rovno91Data(); // 18. m, 1991
TGraphErrors *rovno88Data(); // 18., 25. m, 1988
TGraphErrors *goesgenData(); // 38., 45., et 65. m, 1995
TGraphErrors *bugey3Data(); // 15., 40., et 95. m, 1995
TGraphErrors *ILL1980Data(); // 8.76 m, 1980
TGraphErrors *bugey4Data(); // 15. m, 1994

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *kamlandData(){
// PhysRevLett.108.171803

 

  const char name[]="kamland";
  const double distance[]={180.e3}; // en metre ??
  const double ratioNew[]={0.658}; // 2005
  const double eratiomNew[]={0.044}; // 2005
  const double eratiopNew[]={0.047}; // 2005


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	double eratiom=eratiomNew[i];
	double eratiop=eratiopNew[i];
	printf("%s : %f m  R= %f +- %f (stat) + %f (syst)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *renoData(){
// PhysRevLett.108.171803

 

  const char name[]="RENO";
  const double distance[]={1383}; // en metre
  const double ratioNew[]={0.929}; // arxiv:hep-ex/1312.4111 Proceedings of XVth International Workshop on Neutrino Telescopes (March 2013 at Venice, Italy)
  const double eratiomNew[]={0.006}; //sys
  const double eratiopNew[]={0.007}; //stat


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	// reactor anomaly : normalization parameter a=6% ( nuTurn 2012 http://agenda.infn.it/contributionListDisplay.py?confId=4722 )
	/* ratio/=1.06; */
	double eratiom=eratiomNew[i];
	double eratiop=eratiopNew[i];
	printf("%s : %f m  R= %f +- %f (stat) + %f (syst)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *dayabayData(){
// PhysRevLett.108.171803

 

  const char name[]="dayabay";
  const double distance[]={573.}; // en metre ??
  const double ratioNew[]={0.944}; // valeur de mai 2012 (neutrino 2012) rapport far/near http://neu2012.kek.jp/neu2012/programme.html
  const double eratiomNew[]={0.007};
  const double eratiopNew[]={0.003};


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetLineColor(8);
	 gr->SetMarkerColor(8);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	double eratiom=eratiomNew[i];
	double eratiop=eratiopNew[i];
	printf("%s : %f m  R= %f +- %f (stat) + %f (syst)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *pauloVerdeData(){
// F. Boehm Phys. Rev. D64, 11 (2001).

 

  const char name[]="pauloVerde";
  const double distance[]={820.}; // en metre ??
  const double ratioNew[]={0.975};
  const double eratiomNew[]={0.023};
  const double eratiopNew[]={0.055};


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	double eratiom=eratiomNew[i];
	double eratiop=eratiopNew[i];
	printf("%s : %f m  R= %f +- %f (stat) + %f (syst+Stot)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *choozData(){
// M. Apollonio et al., Phys. Lett. B466, 415 (1999). M. Apollonio et al., Eur. Phys. J. C27, 331-374 (2003).



  const char name[]="chooz";
  const double distance[]={1050.}; // en metre
  const double ratioNew[]={0.961};
  const double eratiomNew[]={0.027};
  const double eratiopNew[]={0.032};


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	double eratiom=eratiomNew[i];
	double eratiop=eratiopNew[i];
	printf("%s : %f m  R= %f +- %f (stat) + %f (syst+Stot)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *SRP96Data(){
// Savannah River Plant (SRP)
// C.L. Jr. Cowan, F. Reines, F.B. Harrison, H.W. Kruse and A.D. McGuire, Science 124 (1956) 103
// Z.D. Greenwood et al., Phys. Rev. D53 11 (1996).

 

  const char name[]="SRP";
  const double distance[]={18.,24.}; // en metre
  const double ratioNew[]={0.953,1.019};
  const double eratiomNew[]={0.006,0.010};
  const double eratiopNew[]={0.0353,0.0377};


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	double eratiom=eratiomNew[i];
	double eratiop=eratiopNew[i];
	printf("%s : %f m  R= %f +- %f (stat) + %f (syst+Stot)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *krasnoyarskData(){
// G.S. Vidyakin et al., JETP. 93 (1987) 424-431.
// G.S. Vidyakin et al., JETP Lett. 59 (1994) 390.

  const char name[]="krasnoyarsk";
  const double distance[]={33.,92.,57.}; // en metre
  const double mesure[]={6.19,6.30,6.26}; // 10^-43 cm2/fission
  const double error[] ={0.36,1.28,0.26};

  const double predNew[] ={6.61,6.61,6.61};
  //const double epredNew[]={0.21,0.21,0.21};
  const double epredNew[]={0.17,0.17,0.17}; // !!! les erreurs utilisees ne sont pas explicitement mises dans l'article. un oublie?

  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=mesure[i]/predNew[i];
	double eratiom=error[i]/predNew[i];
	double eratiop=epredNew[i]*mesure[i]/(predNew[i]*predNew[i]);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	printf("%s : %f m  R= %f +- %f (stat+syst) +- %f (Stot)  (%f total)\n",name,distance[i],ratio,eratiom,eratiop,eratio);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *rovno91Data(){
//  V. Kuvshinnikov et al., JETP 54 N5 (1991) 259.

  const char name[]="rovno91";
  const double distance[]={18.}; // en metre
  const double ratioNew[]={0.940};
  const double eratiomNew[]={0.036};
  const double eratiopNew[]={0.};


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	double eratiom=eratiomNew[i];
	double eratiop=eratiopNew[i];
	printf("%s : %f m  R= %f +- %f (stat+syst+Stot)\n",name,distance[i],ratio,eratiom);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *rovno88Data(){
// A.I. Afonin et al., JETP 94 (1988) 1-17.

  const char name[]="rovno88";
  const double distance[]={18.,18.,18.,25.,18.}; // en metre
  const double ratioNew[]={0.917,0.948,0.972,0.959,0.938};
  const double eratiomNew[]={6.9/100.,6.9/100.,7.8/100.,7.8/100.,7.8/100.};
  const double eratiopNew[]={6.9/100.,6.9/100.,7.2/100.,7.2/100.,7.2/100.};


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	double eratiom=eratiomNew[i]*ratioNew[i];
	double eratiop=eratiopNew[i]*ratioNew[i];
	printf("%s : %f m  R= %f +- %f (stat) +- %f (syst+Stot)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *goesgenData(){
// G. Zacek et al., Phys. Rev. D34 (1986) 2621.

  const char name[]="Goesgen";
  const double distance[]={38.,45.,65.}; // en metre
  const double ratioNew[]={0.966,0.991,0.924};
  const double eratiomNew[]={0.017,0.019,0.033};
  const double eratiopNew[]={0.060,0.062,0.062};


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	double eratiom=eratiomNew[i];
	double eratiop=eratiopNew[i];
	printf("%s : %f m  R= %f +- %f (stat) +- %f (syst+Stot)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *bugey3Data(){
// B. Achkar et al., Nucl. Phys. B434 (1995) 503.

  const char name[]="bugey3";
  const double distance[]={15.,40.,95.}; // en metre
  const double ratioNew[]={0.946,0.952,0.876};
  const double eratiomNew[]={0.004,0.01,0.126};
  const double eratiopNew[]={0.048,0.048,0.048};


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	double eratiom=eratiomNew[i];
	double eratiop=eratiopNew[i];
	printf("%s : %f m  R= %f +- %f (stat) +- %f (syst+Stot)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *ILL1980Data(){
// H. Kwon et al., Phys. Rev. D24 1097 (1981).

  const char name[]="ILL80";
  const double distance[]={8.76}; // en metre
//   const double ratioOld[]={0.832};
//   const double eratiomOld[]={0.029};
//   const double eratiopOld[]={0.0738};

  const double ratioNew[]={0.802};
  const double eratiomNew[]={0.028};
  const double eratiopNew[]={0.071};


  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=ratioNew[i];
	double eratiom=eratiomNew[i];
	double eratiop=eratiopNew[i];
	printf("%s : %f m  R= %f +- %f (stat) +- %f (syst+Stot)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors *bugey4Data(){
// Y. D´clais et al., Phys. Lett. B338 (1994) 383.

  const char name[]="Bugey4";
  const double distance[]={15.}; // en metre
  const double mesure[]={5.752};
  const double error[] ={1.4/100.};

//   const double predOld[] ={5.824};
//   const double epredOld[]={2.7/100.};

  const double predNew[] ={6.102};
  const double epredNew[]={2.7/100.};

  TGraphErrors *gr = new TGraphErrors();
	 gr->SetMarkerColor(1);
	 gr->SetMarkerStyle(20);
	 gr->SetMarkerSize(0.5);
	 gr->SetName(name);

  for(unsigned int i=0;i<sizeof(distance)/sizeof(double);i++){
	double ratio=mesure[i]/predNew[i];
	double eratiom=error[i]*mesure[i]/predNew[i];
	double eratiop=epredNew[i]*predNew[i]*mesure[i]/(predNew[i]*predNew[i]);
	printf("%s : %f m  R= %f +- %f (stat+syst) +- %f (Stot)\n",name,distance[i],ratio,eratiom,eratiop);
	double eratio=sqrt(eratiom*eratiom+eratiop*eratiop);
	gr->SetPoint(gr->GetN(),distance[i],ratio);
	gr->SetPointError(gr->GetN()-1,0.,eratio);
  }
  //printf("%d points\n",gr->GetN());
  return gr;
}

#endif
