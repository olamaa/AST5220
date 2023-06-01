#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  Vector k_array = Utils::linspace(k_min, k_max, n_k);

  Vector k_array_0 = Utils::linspace(k_min_0, k_max_0, n_k_0);
  Vector k_array_1 = Utils::linspace(k_min_1, k_max_1, n_k_1);
  Vector k_array_2 = Utils::linspace(k_min_2, k_max_2, n_k_2);
  Vector k_array_3 = Utils::linspace(k_min_3, k_max_3, n_k_3);
  Vector k_array_4 = Utils::linspace(k_min_4, k_max_4, n_k_4);
  Vector k_array_5 = Utils::linspace(k_min_5, k_max_5, n_k_5);

  // gather the different areas into one single
  all_k_arrays.push_back(k_array_0);
  all_k_arrays.push_back(k_array_1);
  all_k_arrays.push_back(k_array_2);
  all_k_arrays.push_back(k_array_3);
  all_k_arrays.push_back(k_array_4);
  all_k_arrays.push_back(k_array_5);

  //=========================================================================
  // DONE: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  std::cout << "line of sight integration" << std::endl;
  line_of_sight_integration(k_array);

  //=========================================================================
  // DONE: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());

  double dz = 2.*M_PI/10.;
  int n_z = round(k_max*cosmo->eta_of_x(0.)/dz);
  Vector z_array = Utils::linspace(0.,k_max*cosmo->eta_of_x(0.),1e3);
  Vector bessel;
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

  for (int l=0;l < ells.size();l++){
    int ell = ells[l];
    bessel.clear();
    for (int z_index=0; z_index<z_array.size(); z_index++){
      int z = z_array[z_index];
      bessel.push_back(Utils::j_ell(ell, z));
    }
    j_ell_splines[l].create(z_array,bessel,"bessel");
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result;
  Vector ell_result;
  const double eta_0 = cosmo->eta_of_x(0.);
  const double Mpc = Constants.Mpc; 


  Vector current_k_array;
  double dx;
  double area;
  double k;
  double npts;

  //=============================================================================
  // DONE: Implement to solve for the general line of sight integral 
  // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
  // given value of k
  //=============================================================================
  std::cout << "l (of 60)" << "         " << "ell (of 2000)" << "         " << "Percent completed" << std::endl;
  for (int l = 0;l<ells.size();l++){
    double ell = ells[l];
    double variable = l;
    // Where ell is small, we are looking at large scales, and dx should be small 
    std::cout << l << ":                 " << ell << ":                     " << variable/ells.size()*100 << "% " << std::endl;
    if (ell <=6){
      current_k_array = all_k_arrays[0];
      npts = n_k_0;
    }
    if (ell > 6 && ell<=20){
      current_k_array = all_k_arrays[1];
      npts = n_k_1;
    }

    if (ell > 20 && ell<=200){
      current_k_array = all_k_arrays[2];
      npts = n_k_2;
    }
    if (ell > 200 && ell<=400){
      current_k_array = all_k_arrays[3];
      npts = n_k_3;
    }
    if (ell > 400 && ell<=800){
      current_k_array = all_k_arrays[4];
      npts = n_k_4;
    }
    if (ell > 800){
      current_k_array = all_k_arrays[5];
      npts = n_k_5;
    }

    for(int ik = 0; ik < npts; ik++){
      area = 0.0;
      k = current_k_array[ik];
      double bessel_0;
      double bessel_1;
      double ck = Constants.c*k;
      Vector x_array_source;
      x_array_source.clear();
      x_array_source.push_back(-8.);

      for (int i=0; i<1e9;i++){
        bessel_0 = Utils::j_ell(ell, (k*(eta_0-cosmo->eta_of_x(x_array_source[i]))));
        dx = 2.*M_PI*cosmo->Hp_of_x(x_array_source[i])/(10.*10.*10.*ck);
        if (ell <= 10){ // dx should be small for small l's
            dx = 2.*M_PI*cosmo->Hp_of_x(x_array_source[i])/(10.*10.*10.*10.*ck);
        }
        if (ell > 100){ // dx is allowed to be large for larger l's
            dx = 2.*M_PI*cosmo->Hp_of_x(x_array_source[i])/(10.*10.*ck);
        }
        // Using the trapzoidal rule, we get
        if (x_array_source[i]+dx >= 0.){
          x_array_source.push_back(0.);
          bessel_1 = Utils::j_ell(ell, (k*(eta_0-cosmo->eta_of_x(x_array_source[i+1]))));
          dx = 0.- x_array_source[i];
          area += ((source_function(x_array_source[i+1],k)*bessel_1)+(source_function(x_array_source[i],k)*bessel_0))/2.*dx;  
          
          break;
        }
        x_array_source.push_back(x_array_source[i]+dx);
        bessel_1 = Utils::j_ell(ell, (k*(eta_0-cosmo->eta_of_x(x_array_source[i+1]))));
        area += ((source_function(x_array_source[i+1],k)*bessel_1)+(source_function(x_array_source[i],k)*bessel_0))/2.*dx;
      }
      if (x_array_source[x_array_source.size()-1] < 0.){
        std::cout << "x did not reach 0 in LOS" << std::endl;
      }
      ell_result.push_back(area);

    }
    result.push_back(ell_result);
    ell_result.clear();
  }
  std::cout<< "Line of sight integration completed" <<std::endl;
  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// DONE: The line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // DONE: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);
  //std::cout << "Line 231" << std::endl;
  // Spline the result and store it in thetaT_ell_of_k_spline
  for (int i=0;i<nells;i++){
    int ell = ells[i];
    if (ell<=6){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[0], thetaT_ell_of_k[i]); 
    }
    if (ell > 6 && ell <= 20){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[1], thetaT_ell_of_k[i]); 
    }
    if (ell > 20 && ell <= 200){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[2], thetaT_ell_of_k[i]); 
    }
    if (ell > 200 && ell <= 400){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[3], thetaT_ell_of_k[i]); 
    }
    if (ell > 400 && ell <= 800){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[4], thetaT_ell_of_k[i]); 
    }
    if (ell > 800){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[5], thetaT_ell_of_k[i]); 
    }
    std::cout << i << std::endl;
  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    std::vector<Spline> & f_ell_spline){
  std::cout<< "Solving for cell" << std::endl;
  const int nells      = ells.size();
  const double eta_0 = cosmo->eta_of_x(0.);
  double area;
  double dk;
  double dlogk;
  double k_min_cell;
  double k_max_cell;
  Vector result;
  result.clear();
  for (int l = 0; l<nells; l++){
    double ell = ells[l];
    if (ell <=6){
      k_min_cell = all_k_arrays[0][0];
      k_max_cell = all_k_arrays[0][all_k_arrays[0].size()-1];
    }
    if (ell > 6 && ell<=20){
      k_min_cell = all_k_arrays[1][0];
      k_max_cell = all_k_arrays[1][all_k_arrays[1].size()-1];
    }
    if (ell > 20 && ell<=200){
      k_min_cell = all_k_arrays[2][0];
      k_max_cell = all_k_arrays[2][all_k_arrays[2].size()-1];
    }
    if (ell > 200 && ell<=400){
      k_min_cell = all_k_arrays[3][0];
      k_max_cell = all_k_arrays[3][all_k_arrays[3].size()-1];
    }
    if (ell > 400 && ell<=800){
      k_min_cell = all_k_arrays[4][0];
      k_max_cell = all_k_arrays[4][all_k_arrays[4].size()-1];
    }
    if (ell > 800){
      k_min_cell = all_k_arrays[5][0];
      k_max_cell = all_k_arrays[5][all_k_arrays[5].size()-1];
    }
    //=============================================================================
    // DONE: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    //std::cout<<"before dlogk (306)"<<std::endl;
    area = 0.0;
    Vector current_k_array = Utils::linspace(k_min_cell,k_max_cell,1e5);
    arma::vec log_k_array_cell = log(current_k_array);
    for (int index = 0;index < current_k_array.size()-1;index++){
      dlogk = log_k_array_cell[index+1]-log_k_array_cell[index];

      double pps_k1 = primordial_power_spectrum(current_k_array[index+1])*pow(f_ell_spline[l](current_k_array[index+1]),2);
      double pps_k = primordial_power_spectrum(current_k_array[index])*pow(f_ell_spline[l](current_k_array[index]),2);
      area += (pps_k1 + pps_k)/2. * dlogk;
    }
    result.push_back(4.*M_PI*area);
  }
  //std::cout<<"after dlogk (310)"<<std::endl;
  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k) const{
  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================
  double delta_m = Constants.c*Constants.c*k*k*pert->get_Phi(x,k)*2./(3.*(cosmo->get_OmegaB(0.)+cosmo->get_OmegaCDM(0.))*cosmo->get_H0()*cosmo->get_H0())*exp(x);
  double pofk = pow(abs(delta_m),2)*primordial_power_spectrum(k)*2.*M_PI*M_PI/(k*k*k);

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorell = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

// outputs theta
void PowerSpectrum::output_thetak(const std::string filename) const{
  const int    n_pts =  10000;
  const double Mpc = Constants.Mpc;
  const double c_H0 = Constants.c/cosmo->get_H0();
  const double norm = 1e6*cosmo->get_H0()/Constants.c;
  std::cout << "norm: " << norm << std::endl;
  arma::vec k_array = arma::logspace(log10(k_min),log10(k_max),n_pts);
  std::ofstream fp(filename.c_str());
  const int size_of_theta_ell = thetaT_ell_of_k_spline.size();
  std::cout << size_of_theta_ell << std::endl;
  // for the different regimes of ell
  auto print_data = [&] (const double k) {
    fp <<    k*c_H0                                     << " ";
    fp <<    k*norm                                     << " "; // 10**-47..? 
    for(int i=0; i < size_of_theta_ell; i++){
      fp << thetaT_ell_of_k_spline[i](k)                  << " ";   // ell
    };
    fp <<"\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}


void PowerSpectrum::output_matter(const std::string filename) const{
  const int    n_pts =  10000;
  const double Mpc = Constants.Mpc;
  const double conversion_P = pow(cosmo->get_h(),3)/pow(Mpc,3); // Conversion to Mpc
  const double conversion = Mpc/cosmo->get_h();
  Vector k_array = Utils::linspace(k_min, k_max, n_pts);
  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double k) {
    fp <<    k*conversion                                     << " ";
    fp << get_matter_power_spectrum(0.,k)*conversion_P        << " ";
    fp <<"\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}


void PowerSpectrum::output_Source(const std::string filename) const{
  const double x_min_ = -8.;
  const double x_max_ =  0.;
  const int    n_pts =  10000;
  double k = 340.*cosmo->get_H0()/Constants.c;
  const double Mpc = Constants.Mpc;
  Vector x_array = Utils::linspace(x_min_,x_max_,n_pts);
  double eta_0 = cosmo->eta_of_x(0.);
  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp <<    x                                                                           << " ";
    fp << pert->get_Source_T(x,k)*Utils::j_ell(100,(k*(eta_0-cosmo->eta_of_x(x))))* 1e3  << " "; // j_ell(100, k...) where 100 is ell to compare with Callin (2006)
    fp << pert->get_Source_T(x,k)                                                        << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}