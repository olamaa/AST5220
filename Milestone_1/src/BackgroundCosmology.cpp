 #include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{
  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ... // DONE
  //=============================================================================
  H0 = Constants.H0_over_h*h;
  OmegaR = 2. * M_PI*M_PI/30. * pow(Constants.k_b*TCMB, 4)/(pow(Constants.hbar, 3)*pow(Constants.c, 5)) * 8.*M_PI*Constants.G/(3.*pow(H0, 2));
  OmegaNu = Neff * (7./8.) * pow(4./11., 4./3.) * OmegaR;
  OmegaLambda = 1.-(OmegaK + OmegaB + OmegaCDM + OmegaR + OmegaNu);
  //std::cout << OmegaLambda << std::endl;
}


//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines // DONE
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  int npts = 1e4;

  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    detadx[0] = Constants.c/Hp_of_x(x);

    return GSL_SUCCESS;
  };
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    //=============================================================================
    // TODO: Set the rhs of the dtdx ODE
    //=============================================================================
    dtdx[0] = 1/H_of_x(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make // DONE
  // the spline eta_of_x_spline 
  //=============================================================================
  Vector eta0{0.};
  ODESolver ode;
  ode.solve(detadx, x_array, eta0);
  auto eta_array = ode.get_data_by_component(0); // try c/(Hp(x=x_start)), else: 0
  eta_of_x_spline.create(x_array,eta_array,"eta_of_x");



  Vector t0{1/(2*H_of_x(x_start))}; // Initial condition
  ode.solve(dtdx, x_array, t0);
  auto t_array = ode.get_data_by_component(0);
  t_of_x_spline.create(x_array,t_array,"t_of_x");

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{

  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================
  double exponentials = (OmegaB+OmegaCDM)*pow(exp(x), -3) + (OmegaR + OmegaNu)*pow(exp(x), -4) + OmegaK*pow(exp(x), -2) + OmegaLambda;
  double H = H0*sqrt(exponentials);
  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{

  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================
  double Hp = exp(x)*H_of_x(x);
  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================

  // YOU NEED TO FIX THIS 

  double exponential = (OmegaB+OmegaCDM)*pow(exp(x),-1)+ (OmegaR+OmegaNu)*pow(exp(x),-2)+OmegaK+OmegaLambda*pow(exp(x),2);
  double d_exponential = -(OmegaB+OmegaCDM)*pow(exp(x),-1) -2.* (OmegaR+OmegaNu)*pow(exp(x),-2.)+2.*OmegaLambda*pow(exp(x),2);
  return H0*1./2.*exponential*d_exponential;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================

  // YOU NEED TO FIX THIS 
  
  double exponential = (OmegaB+OmegaCDM)*pow(exp(x),-1.) + (OmegaR+OmegaNu)*pow(exp(x),-2.)+OmegaK + OmegaLambda*pow(exp(x),2.);
  double d_exponential = -(OmegaB+OmegaCDM)*pow(exp(x),-1.) - 2.* (OmegaR+OmegaNu)*pow(exp(x),-2.) + 2.*OmegaLambda*pow(exp(x),2.);
  double dd_exponential = (OmegaB+OmegaCDM)*pow(exp(x),-1.) + 4.* (OmegaR+OmegaNu)*pow(exp(x),-2.) + 4.*OmegaLambda*pow(exp(x),2.);

  return H0*1./2.*(d_exponential*d_exponential+exponential*dd_exponential);
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================
  return OmegaB*H0*H0/(pow(exp(x), 3) * H_of_x(x)*H_of_x(x));
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================
  return OmegaR*H0*H0/(pow(exp(x), 4.)*H_of_x(x)*H_of_x(x));
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================
  return OmegaNu*H0*H0/(pow(exp(x), 4)*H_of_x(x)*H_of_x(x));
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================
  return OmegaCDM*H0*H0/(pow(exp(x), 3.)*H_of_x(x)*H_of_x(x));
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================

  return OmegaLambda*H0*H0/(H_of_x(x)*H_of_x(x));
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================
  return OmegaK*H0*H0/(pow(exp(x), 2)*H_of_x(x)*H_of_x(x));
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================
  if (OmegaK < 0.){
  double r = get_comoving_distance_of_x(x)*sin(sqrt(-OmegaK)*H0*get_comoving_distance_of_x(x)/Constants.c) / (sqrt(-OmegaK)*H0*get_comoving_distance_of_x(x)/Constants.c);
    return exp(-x)*r;
  }
  else if (OmegaK == 0.){
  return exp(-x)*get_comoving_distance_of_x(x);
  }
  else if (OmegaK > 0.){
  double r = get_comoving_distance_of_x(x)*sinh(sqrt(OmegaK)*H0*get_comoving_distance_of_x(x)/Constants.c) / (sqrt(OmegaK)*H0*get_comoving_distance_of_x(x)/Constants.c);
  return exp(-x)*r;
  }
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement... // DONE
  //=============================================================================
  return eta_of_x(0) - eta_of_x(x);
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB                           << "    \n";
  std::cout << "OmegaCDM:    " << OmegaCDM                         << "    \n";
  std::cout << "OmegaLambda: " << OmegaLambda                      << "    \n";
  std::cout << "OmegaK:      " << OmegaK                           << "    \n";
  std::cout << "OmegaNu:     " << OmegaNu                          << "    \n";
  std::cout << "OmegaR:      " << OmegaR                           << "    \n";
  std::cout << "Neff:        " << Neff                             << "    \n";
  std::cout << "h:           " << h                                << "    \n";
  std::cout << "TCMB:        " << TCMB                             << "    \n";
  std::cout << "time:        " << t_of_x(0)/(1e9*365*24*60*60)  << "Gyr \n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -20.0;
  const double x_max =  5.0;
  const int    n_pts =  10000;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                                    << " ";
    fp << eta_of_x(x)                          << " ";
    fp << Hp_of_x(x)                           << " ";
    fp << dHpdx_of_x(x)                        << " ";
    fp << get_OmegaB(x)                        << " ";
    fp << get_OmegaCDM(x)                      << " ";
    fp << get_OmegaLambda(x)                   << " ";
    fp << get_OmegaR(x)                        << " ";
    fp << get_OmegaNu(x)                       << " ";
    fp << get_OmegaK(x)                        << " ";
    fp << get_luminosity_distance_of_x(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

