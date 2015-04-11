// ************** Gauss-Legendre Integrator ************** //
void gauleg(const double x1, const double x2, int n, double *x, double *w);

// ************** Compute the Neutrinos flux in the Beam from ILL reactor ************** //
double neutrino_flux(double E_ve);

// ************** Compute antineutrino cross-section  ************** //
double neutrino_sigma(double E_ve);
