#include <math.h>

// ************** Gauss-Legendre Integrator ************** //
void gauleg(const double x1, const double x2, int n, double *x, double *w)
{
    const double EPS=1.0e-14;
    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;

    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for (i=0;i<m;i++) {
        z=cos(3.141592654*(i+0.75)/(n+0.5));
        do {
            p1=1.0;
            p2=0.0;
            for (j=0;j<n;j++) {
                p3=p2;
                p2=p1;
                p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        } while (fabs(z-z1) > EPS);
        x[i]=xm-xl*z;
        x[n-1-i]=xm+xl*z;
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n-1-i]=w[i];
    }
} 

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
