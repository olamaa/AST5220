#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();

  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
  sound_horizon();

}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on // DONE
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr = Utils::linspace(0., 1., npts_rec_arrays);
  Vector ne_arr = Utils::linspace(0., 1., npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;

  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so // DONE
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation // DONE
      //=============================================================================
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode // DONE
      //==============================================================

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        dXedx[0] = Xe_saha_limit;
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result //DONE
      //=============================================================================


      Vector Xe_init{Xe_current};
      const double OmegaB0      = cosmo->get_OmegaB(0.0);
      const double G           = Constants.G;

      const double rho_c0 = 3.*pow(cosmo->get_H0(),2.)/(8*M_PI*G);
      const double m_H         = Constants.m_H;

      Vector x_ = {x_array.begin()+i, x_array.end()};   
      peebles_Xe_ode.solve(dXedx, x_, Xe_init);
      auto Xe_p = peebles_Xe_ode.get_data_by_component(0);
      for (int j=i; j < npts_rec_arrays;j++){
        Xe_arr[j] = Xe_p[j-i];
        double nb = OmegaB0*rho_c0/(m_H*exp(3.*x_array[j]));
        ne_arr[j] = Xe_p[j-i]*nb;

      }
      break;
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working // DONE
  //=============================================================================



  Vector log_ne = Utils::linspace(0.,1,npts_rec_arrays); 
  for (int i=0;i<npts_rec_arrays;i++){
    log_ne[i] = log(ne_arr[i]);
  }
  
  Xe_of_x_spline.create(x_array,Xe_arr,"Xe");
  log_ne_of_x_spline.create(x_array,log_ne,"ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double OmegaB0      = cosmo->get_OmegaB(0.0);
  double TCMB               = cosmo->get_TCMB(x);

  // Electron fraction and number density
  double Xe;
  double ne;
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================


  double Tb = TCMB;
  const double rho_c0 = 3*pow(cosmo->get_H0(), 2)/(8*M_PI*G);
  double nb = OmegaB0*rho_c0/(m_H*pow(a, 3));
  double nH = (1 - Yp)*nb; // nH = nb
  double b = 1./nb * pow(m_e*Tb*k_b/(2.*M_PI*hbar*hbar), 3./2.)*exp(-epsilon_0/(Tb*k_b));
  Xe =  2./(sqrt(1.+4./b) + 1.);
  ne = Xe*nH;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  const double rho_c0       = 3*pow(cosmo->get_H0(), 2)/(8*M_PI*G);
  const double OmegaB0      = cosmo->get_OmegaB(0.0);
  double nb                 = OmegaB0*rho_c0/(m_H*pow(a, 3));
  double TCMB               = cosmo->get_TCMB(x);
  double H                  = cosmo->H_of_x(x);
  //=============================================================================
  // TODO: Write the expression for dXedx // DONE
  //=============================================================================

  double Tb = TCMB;
  double phi_2       = 0.448*log(epsilon_0/(TCMB*k_b));                                       
  double alpha_2     = c*8./sqrt(3.*M_PI)*sigma_T*sqrt(epsilon_0/(TCMB*k_b))*phi_2;       
  
  double beta        = alpha_2*pow(m_e*TCMB*k_b/(2.*M_PI*hbar*hbar), 3./2.)*exp(-epsilon_0/(TCMB*k_b));
  double beta_2      = alpha_2*pow(m_e*TCMB*k_b/(2.*M_PI*hbar*hbar), 3./2.)*exp(-1./4.*epsilon_0/(TCMB*k_b));                                      

  double nH          = (1. - Yp)*nb;                                                        

  double n_1s        = (1. - X_e)*nH;                                                       
  double lambda_a    = H * pow(3.*epsilon_0/(hbar*c), 3.)/(pow(8.*M_PI, 2.)*n_1s);                    


  //std::cout << X_e << std::endl;

  double Cr          = (lambda_2s1s + lambda_a)/(lambda_2s1s + lambda_a + beta_2);
  double rhs         = Cr/H*(beta*(1.-X_e) - nH*alpha_2*X_e*X_e);                         
  //std::cout <<  << std::endl;
  dXedx[0] = rhs;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result // DONE
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 10000;
  Vector x_array = Utils::linspace(x_start, x_today, npts);
  Vector x_array_r = Utils::linspace(x_today, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    //=============================================================================
    // TODO: Write the expression for dtaudx // DONE
    //=============================================================================
    const double c = Constants.c;
    const double sigma_T = Constants.sigma_T;
    const double H = cosmo->H_of_x(x);
    const double ne = ne_of_x(x);

    dtaudx[0] = -c*ne*sigma_T/H;
    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines // DONE
  //=============================================================================

  ODESolver tau_of_x_ode;
  Vector tau_init = {0.0};
  tau_of_x_ode.solve(dtaudx,x_array_r,tau_init);
  auto tau_of_x_solution = tau_of_x_ode.get_data_by_component(0);
  Vector tau_of_x_linspace = Utils::linspace(0., 1., npts);
  for (int i=0;i<npts;i++){
    tau_of_x_linspace[i] = tau_of_x_solution[npts - 1 - i];
  }
  tau_of_x_spline.create(x_array,tau_of_x_linspace,"tau");



  //=============================================================================
  // TODO: Compute visibility functions and spline everything // DONE
  //=============================================================================
  Vector g_tilde_linspace = Utils::linspace(0., 1., npts);
  for (int i=0;i<npts;i++){
    g_tilde_linspace[i] = g_tilde(x_array[i]);
  }
  g_tilde_of_x_spline.create(x_array, g_tilde_linspace, "g");
  Utils::EndTiming("opticaldepth");
}

double RecombinationHistory::g_tilde(double x) const{
  return -dtaudx_of_x(x)*exp(-tau_of_x(x));
}


void RecombinationHistory::sound_horizon(){
  Vector x_array = Utils::linspace(x_early, x_today, 2*npts_rec_arrays);
  const double c = Constants.c;
  const double sigma_T = Constants.sigma_T;
  const double OmegaR0 = cosmo->get_OmegaR(0.0);
  const double OmegaB0 = cosmo->get_OmegaB(0.0);

  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
    const double a = exp(x);
    const double Hp = cosmo->Hp_of_x(x);
    const double R = 4.*OmegaR0/(3.*OmegaB0*a);
    const double cs = c*sqrt(R/(3.*(1. + R)));
    dsdx[0] = cs/Hp;
    return GSL_SUCCESS;
  };
  const double R_initial = 4.*OmegaR0/(3*OmegaB0*exp(x_start));
  const double cs_initial = c*sqrt(R_initial/(3.*(1. + R_initial)));
  const double Hp_initial = cosmo->Hp_of_x(x_start);

  ODESolver sound_horizon_ODE;
  Vector sound_horizon_initial = {cs_initial/Hp_initial};
  sound_horizon_ODE.solve(dsdx, x_array, sound_horizon_initial);
  auto sound_horizon_solution = sound_horizon_ODE.get_data_by_component(0);
  sound_horizon_spline.create(x_array, sound_horizon_solution, "SH");
}

//====================================================
// Get methods // DONE
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}
double RecombinationHistory::get_sound_horizon() const{
  return sound_horizon_spline(x_decoupling);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this // DONE
  //=============================================================================

  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement //DONE
  //=============================================================================

  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement // DONE
  //=============================================================================

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement // DONE
  //=============================================================================

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement // DONE
  //=============================================================================

  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================

  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:                              " << Yp                                                         << "\n";
  std::cout << "Xe at decoupling:                " << std::setprecision(15) << Xe_of_x(x_decoupling)             << "\n";
  std::cout << "Sound Horizon [Gpc]:             " << get_sound_horizon()/(Constants.Mpc*1e3)                    << "\n";
  std::cout << "Sound Horizon: conformal ratio:  " << get_sound_horizon()/cosmo->eta_of_x(x_decoupling)          << "\n";
  std::cout << "\n";
  std::cout << "x_decoupling:                    " << std::setprecision(15) << x_decoupling                      << "\n";
  std::cout << "tau_decoupling:                  " << std::setprecision(15) << tau_of_x(x_decoupling)            << "\n"; 
  std::cout << "z_decoupling:                    " << std::setprecision(15) << exp(-x_decoupling)-1              << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_today;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp                          << x                     << " ";
    fp << std::setprecision(15) << Xe_of_x(x)            << " "; 
    fp                          << ne_of_x(x)            << " "; 
    fp                          << tau_of_x(x)           << " "; 
    fp                          << dtaudx_of_x(x)        << " "; 
    fp                          << ddtauddx_of_x(x)      << " "; 
    fp                          << g_tilde_of_x(x)       << " "; 
    fp                          << dgdx_tilde_of_x(x)    << " "; 
    fp                          << ddgddx_tilde_of_x(x)  << " "; 
    fp                          << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

