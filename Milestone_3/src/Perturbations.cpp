#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo,  
    RecombinationHistory *rec) : 
  cosmo(cosmo),  
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){                                                                                       

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  //compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){                                                                                     
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // DONE: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  arma::vec k_array_arma = arma::logspace(log10(k_min),  log10(k_max),  n_k);
  Vector k_array;
  for (int i=0;i<k_array_arma.size();i++){
      k_array.push_back(k_array_arma[i]);
  }
  //std::cout << "Line 45" << std::endl;
  // Loop over all wavenumbers

  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
    if(ik == n_k-1) std::cout << std::endl;
    }
  std::cout << "Progress bar worked" << std::endl;

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);

    std::cout << "k and tight coupling time sat" << std::endl;

    //===================================================================
    // DONE: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    int tight_index;
    for (int i =0;i < x_array.size();i++){
      if (x_array[i] > x_end_tight){
        tight_index = i;
        break;  
      }
      else{}
    }

    //std::cout << "Line 81" << std::endl;

    Vector tight_coupling_x_array;
    tight_coupling_x_array.clear();
    for (int i=0;i<tight_index;i++){
        tight_coupling_x_array.push_back(x_array[i]);
    }

    //std::cout << "Line 89" << std::endl;

    Vector from_and_after_tc_x_array;
    from_and_after_tc_x_array.clear();

    //std::cout << "Line 94" << std::endl;

    for (int i=0;i < n_x-tight_index+1;i++){
      // Remove last element to avoid overlapping

      //std::cout << "Line 99" << std::endl;

      from_and_after_tc_x_array.push_back(x_array[tight_index-1+i]);

      //std::cout << "Line 103" << std::endl;

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start,  k);                       

    //std::cout << "Line 108" << std::endl;

    // DONE: The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x,  const double *y,  double *dydx){
      return rhs_tight_coupling_ode(x,  k,  y,  dydx);
    };
    //std::cout << "Line 114" << std::endl;

    // DONE: Integrate from x_start -> x_end_tight
    ODESolver tight_coupling_ODE;
    std::vector<Vector> tight_coupling_data;
    //Vector tight_coupling_data;

    //std::cout << "Line 121" << std::endl;

    if (tight_coupling_x_array.size() >= 2){
        tight_coupling_ODE.solve(dydx_tight_coupling, tight_coupling_x_array, y_tight_coupling_ini);
        tight_coupling_data = tight_coupling_ODE.get_data();
    }
    else{
      tight_coupling_data.push_back(y_tight_coupling_ini);
    }

    //std::cout << "Line 131" << std::endl;

    // We only use 7 parameters in tight coupling
    Vector y_tight_coupling_end(7);
    for (int i=0; i < y_tight_coupling_end.size(); i++){
      y_tight_coupling_end[i] = tight_coupling_data.back()[i];
    }
    //std::cout << "Line 138" << std::endl;

    //====i===============================================================
    // DONE: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
     auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling_end,  x_end_tight,  k);

    //std::cout << "Line 150" << std::endl;

    // The full ODE system
    ODEFunction dydx_full = [&](double x,  const double *y,  double *dydx){
      return rhs_full_ode(x,  k,  y,  dydx);
    };

    //std::cout << "Line 157" << std::endl;

    ODESolver full_ODE;
    std::vector<Vector> after_tc_data;

    //std::cout << "Line 162" << std::endl;

    if (from_and_after_tc_x_array.size() >= 2){
      full_ODE.solve(dydx_full, from_and_after_tc_x_array, y_full_ini);
      after_tc_data = full_ODE.get_data();

      //std::cout << "Line 167" << std::endl;

    }
    else{

      //std::cout << "Line 173" << std::endl;

      if (from_and_after_tc_x_array.size() == 1){

        //std::cout << "Line 177" << std::endl;
        //std::cout << tight_coupling_data.back() << std::endl;

        //after_tc_data[0] = {0.};

        //after_tc_data[0] = tight_coupling_data.back();                                                  // HERE ERROR IS HERE !!!!

        std::cout << tight_coupling_data.back().size() << std::endl;

        
        after_tc_data.push_back(tight_coupling_data.back());                    // HERE ERROR HERE HIHI
    

        //std::cout << after_tc_data << std::endl;
        //std::cout << after_tc_data.size() << std::endl;

      }
      else{

        //std::cout << "Line 186" << std::endl;

      }

    //std::cout << "Line 189" << std::endl;

    }

    std::cout << after_tc_data.size() << " " << std::endl;
    //std::cout << "Line 191" << std::endl;

    tight_coupling_data.pop_back();

    //std::cout << "Line 195" << std::endl;


    // DONE: Integrate from x_end_tight -> x_end
    for (int i = 0; i < tight_coupling_data.size() + after_tc_data.size(); i++){
        if (i<tight_coupling_data.size()){
          delta_cdm_vector.push_back(tight_coupling_data[i][0]);
          v_cdm_vector.push_back(tight_coupling_data[i][1]);
          delta_b_vector.push_back(tight_coupling_data[i][2]);
          v_b_vector.push_back(tight_coupling_data[i][3]);
          Phi_vector.push_back(tight_coupling_data[i][4]);
          
          Theta0_vector.push_back(tight_coupling_data[i][5]);
          Theta1_vector.push_back(tight_coupling_data[i][6]);
          //std::cout<< i << std::endl;
          
          //Higher order thetas are zero during tight coupling (approx)
          Theta2_vector.push_back(0.);
          Theta3_vector.push_back(0.);
          Theta4_vector.push_back(0.);
          Theta5_vector.push_back(0.);
          Theta6_vector.push_back(0.);
          Theta7_vector.push_back(0.);
        }
        else{
          delta_cdm_vector.push_back(after_tc_data[i-tight_coupling_data.size()][0]);
          v_cdm_vector.push_back(after_tc_data[i-tight_coupling_data.size()][1]);
          delta_b_vector.push_back(after_tc_data[i-tight_coupling_data.size()][2]);
          v_b_vector.push_back(after_tc_data[i-tight_coupling_data.size()][3]);
          Phi_vector.push_back(after_tc_data[i-tight_coupling_data.size()][4]);

          Theta0_vector.push_back(after_tc_data[i-tight_coupling_data.size()][5]);
          Theta1_vector.push_back(after_tc_data[i-tight_coupling_data.size()][6]);

          Theta2_vector.push_back(after_tc_data[i-tight_coupling_data.size()][7]);
          Theta3_vector.push_back(after_tc_data[i-tight_coupling_data.size()][8]);
          Theta4_vector.push_back(after_tc_data[i-tight_coupling_data.size()][9]);
          Theta5_vector.push_back(after_tc_data[i-tight_coupling_data.size()][10]);
          Theta6_vector.push_back(after_tc_data[i-tight_coupling_data.size()][11]);
          Theta7_vector.push_back(after_tc_data[i-tight_coupling_data.size()][12]);
        }
        ////std::cout << "Line 245" << std::endl;
    }

    // Storing Psi
    
    // Using the initial condition
    Psi_vector.push_back(-2./3.);
    const double OmegaR0 = cosmo->get_OmegaR(0.);
    const double Omegab0 = cosmo->get_OmegaB(0.);
    const double OmegaCDM0 = cosmo->get_OmegaCDM(0.);
    const double H0 = cosmo->get_H0();
    ////std::cout << "Line 256" << std::endl;
    for (int i = 1; i < x_array.size(); i++){
      double a = exp(x_array[i]);
      Psi_vector.push_back(-Phi_vector[i+n_x * ik]-12.*H0*H0/(Constants.c*Constants.c*k*k*a*a)*OmegaR0*Theta2_vector[i + n_x*ik]);
    }
    ////std::cout << "Line 261" << std::endl;



    //===================================================================
    // DONE: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x, k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix,  ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array,  k_array,  f_array);
    // We can now use the spline as f_spline(x,  k)
    //
    // NB: If you use Theta_spline then you have to allocate it first, 
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================


  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // DONE: Make all splines needed: Theta0, Theta1, Theta2, Phi, Psi, ...
  //=============================================================================
  //std::cout << "Line 293" << std::endl;
  std::cout << x_array.size() << " " << k_array.size() << " " << delta_cdm_vector.size() << std::endl;

  delta_cdm_spline.create(x_array,  k_array,  delta_cdm_vector, "delta_cdm_spline");

  //std::cout << "Line 298" << std::endl;

  delta_b_spline.create(x_array, k_array,  delta_b_vector, "delta_b_spline");
  v_cdm_spline.create(x_array, k_array, v_cdm_vector, "v_cdm_spline");
  v_b_spline.create(x_array, k_array, v_b_vector, "v_b_spline");
  Phi_spline.create(x_array, k_array, Phi_vector, "Phi_spline");
  Psi_spline.create(x_array, k_array, Psi_vector, "Psi_spline");

  //std::cout << "Line 302" << std::endl;

  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  Theta_spline[0].create(x_array, k_array, Theta0_vector, "Theta0_spline");
  Theta_spline[1].create(x_array, k_array, Theta1_vector, "Theta1_spline");
  Theta_spline[2].create(x_array, k_array, Theta2_vector, "Theta2_spline");
  Theta_spline[3].create(x_array, k_array, Theta3_vector, "Theta3_spline");
  Theta_spline[4].create(x_array, k_array, Theta4_vector, "Theta4_spline");
  Theta_spline[5].create(x_array, k_array, Theta5_vector, "Theta5_spline");
  Theta_spline[6].create(x_array, k_array, Theta6_vector, "Theta6_spline");
  Theta_spline[7].create(x_array, k_array, Theta7_vector, "Theta7_spline");

  //std::cout << "Line 314" << std::endl;
  }
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x,  const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];    // 1
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];      // 1
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];        // -ck/2H psi
  double &v_b          =  y_tc[Constants.ind_vb_tc];          // -ck/2H psi
  double &Phi          =  y_tc[Constants.ind_Phi_tc];         // 2/3


  //double *Theta        = &y_tc[Constants.ind_start_theta_tc]; // 3
  //double *Nu           = &y_tc[Constants.ind_start_nu_tc];    // 0

  double &Theta0 = y_tc[Constants.ind_start_theta_tc];
  double &Theta1 = y_tc[Constants.ind_start_theta_tc + 1];

  // THE INDEXES FOR THETA AND NU ARE EQUAL

  //=============================================================================
  // DONE: Set the initial conditions in the tight coupling regime
  //=============================================================================
  const double Psi  = -2./3.;
  Phi               = -Psi;
  delta_cdm         = -3./2.*Psi;
  delta_b           = delta_cdm;
  v_cdm             = - Constants.c*k/(2*cosmo->Hp_of_x(x))*Psi;
  v_b               = v_cdm;
  Theta0            = -1./2.*Psi;
  Theta1            = Constants.c*k/(6*cosmo->Hp_of_x(x))*Psi;

  return y_tc;
}
//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc,  
    const double x,  
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double &Theta0_tc        = y_tc[Constants.ind_start_theta_tc];
  const double &Theta1_tc        = y_tc[Constants.ind_start_theta_tc+1];
  //const double Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double &Theta0           = y[Constants.ind_start_theta_tc];
  double &Theta1           = y[Constants.ind_start_theta_tc + 1];
  double &Theta2           = y[Constants.ind_start_theta_tc + 2];
  double &Theta3           = y[Constants.ind_start_theta_tc + 3];
  double &Theta4           = y[Constants.ind_start_theta_tc + 4];
  double &Theta5           = y[Constants.ind_start_theta_tc + 5];
  double &Theta6           = y[Constants.ind_start_theta_tc + 6];
  double &Theta7           = y[Constants.ind_start_theta_tc + 7];
  //double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  //double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // DONE: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  double a = exp(x);
  const double OmegaR0 = cosmo->get_OmegaR(0.0);
  const double Omegab0 = cosmo->get_OmegaB(0.0);
  const double OmegaCDM0 = cosmo->get_OmegaCDM(0.0);
  const double H0 = cosmo->get_H0();
  double dtau_dx = rec->dtaudx_of_x(x);
  double ddtau_ddx = rec->ddtauddx_of_x(x);
  double Hp = cosmo->Hp_of_x(x); 
  double dHp_dx = cosmo->dHpdx_of_x(x);
  double ck_Hp = Constants.c*k/Hp;

  // SET: Scalar quantities (Gravitational potental,  baryons and CDM)
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;
  Phi = Phi_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta0 = Theta0_tc;
  Theta1 = Theta1_tc;
  Theta2 = -20./45.*ck_Hp/dtau_dx*Theta1_tc;

  for (int i = 3; i < n_ell_theta;i++){
    y[Constants.ind_start_theta_tc+i] = -(i/(2.*i+1.))*ck_Hp/dtau_dx*y[Constants.ind_start_theta_tc+i-1];
  }
  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k){
    for (double time_index = 0; time_index < x_array.size(); time_index++){
      if (abs(rec->dtaudx_of_x(x_array[time_index])) > 10. && abs(rec->dtaudx_of_x(x_array[time_index])) > 10.*Constants.c*k/cosmo->Hp_of_x(x_array[time_index])){}
      else{
        if (x_array[time_index] <= -8.3){
          return x_array[time_index];
        }
        else{
          return -8.3;
        }
      }
    } 
    //=============================================================================
  // DONE: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  std::cout << "Invalid value: Time for TC is incorrect." << std::endl;
  return 0;
}


//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // DONE: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix, ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // DONE: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array,  k_array,  ST_array,  "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array,  k_array,  SE_array,  "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x,  double k,  const double *y,  double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double &Theta0           = y[Constants.ind_start_theta_tc];
  const double &Theta1           = y[Constants.ind_start_theta_tc + 1];
  //const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double &dTheta0dx        =  dydx[Constants.ind_start_theta_tc];
  double &dTheta1dx        =  dydx[Constants.ind_start_theta_tc + 1];
  //double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // DONE: fill in the expressions for all the derivatives
  //=============================================================================
  double a = exp(x);
  const double OmegaR0 = cosmo->get_OmegaR(0.0);
  const double Omegab0 = cosmo->get_OmegaB(0.0);
  const double OmegaCDM0 = cosmo->get_OmegaCDM(0.0);
  const double H0 = cosmo->get_H0();
  double dtau_dx = rec->dtaudx_of_x(x);
  double ddtau_ddx = rec->ddtauddx_of_x(x);
  double Hp = cosmo->Hp_of_x(x); 
  double dHp_dx = cosmo->dHpdx_of_x(x);
  double ck_Hp = Constants.c*k/Hp;
  double R = 4.*OmegaR0/(3.*Omegab0*a);
  double Theta2 = -20./45.*ck_Hp/dtau_dx*Theta1;
  double Psi = -Phi-12.*H0*H0/(Constants.c*Constants.c*k*k*a*a)*OmegaR0*Theta2;
  

  dPhidx = Psi-ck_Hp*ck_Hp*Phi/3.+ H0*H0/(2.*Hp*Hp)*(OmegaCDM0/a*delta_cdm+Omegab0/a*delta_b+4.*OmegaR0/(a*a)*Theta0);
  dTheta0dx = -ck_Hp*Theta1-dPhidx;
  double q = (-((1.-R)*dtau_dx+(1.+ R)*ddtau_ddx)*(3.*Theta1+v_b)-ck_Hp*Psi+(1.-dHp_dx/Hp)*ck_Hp*(-Theta0+2.*Theta2)-ck_Hp*dTheta0dx)
             /((1. + R)*dtau_dx+dHp_dx/Hp-1.);

  // SET: Scalar quantities (Phi,  delta,  v,  ...)
  dv_bdx = -v_b-ck_Hp*Psi+R*(1./(1.+R))*(q + v_b + ck_Hp*(-Theta0+2.*Theta2));
  ddelta_cdmdx = ck_Hp*v_cdm-3.*dPhidx;
  dv_cdmdx = -v_cdm-ck_Hp*Psi;
  ddelta_bdx = ck_Hp*v_b-3.*dPhidx;

  // SET: Photon multipoles (Theta_ell)
  dTheta1dx = 1./3.*(q-dv_bdx);

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x,  double k,  const double *y,  double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const int n_ell_tot_full      = Constants.n_ell_tot_full;    
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double &Theta0          =  y[Constants.ind_start_theta_tc];
  const double &Theta1          =  y[Constants.ind_start_theta_tc + 1];
  const double &Theta2          =  y[Constants.ind_start_theta_tc + 2];
  const double &Theta3          =  y[Constants.ind_start_theta_tc + 3];
  const double &Theta4          =  y[Constants.ind_start_theta_tc + 4];
  const double &Theta5          =  y[Constants.ind_start_theta_tc + 5];
  const double &Theta6          =  y[Constants.ind_start_theta_tc + 6];
  const double &Theta7          =  y[Constants.ind_start_theta_tc + 7];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double &dTheta0_dx      = dydx[Constants.ind_start_theta_tc];
  double &dTheta1_dx      = dydx[Constants.ind_start_theta_tc + 1];
  double &dTheta2_dx      = dydx[Constants.ind_start_theta_tc + 2];
  double &dTheta3_dx      = dydx[Constants.ind_start_theta_tc + 3];
  double &dTheta4_dx      = dydx[Constants.ind_start_theta_tc + 4];
  double &dTheta5_dx      = dydx[Constants.ind_start_theta_tc + 5];
  double &dTheta6_dx      = dydx[Constants.ind_start_theta_tc + 6];
  double &dTheta7_dx      = dydx[Constants.ind_start_theta_tc + 7];

  // Cosmological parameters and variables
  double a = exp(x);
  const double OmegaR0 = cosmo->get_OmegaR(0.);
  const double Omegab0 = cosmo->get_OmegaB(0.);
  const double OmegaCDM0 = cosmo->get_OmegaCDM(0.);
  const double H0 = cosmo->get_H0();
  double dtau_dx = rec->dtaudx_of_x(x);
  double ddtau_ddx = rec->ddtauddx_of_x(x);
  double Hp = cosmo->Hp_of_x(x); 
  double dHp_dx = cosmo->dHpdx_of_x(x);
  double eta_x = cosmo->eta_of_x(x);
  double ck_Hp = Constants.c*k/Hp;
  double R = 4.*OmegaR0/(3.*Omegab0*a);
  double Psi = -Phi-12.*H0*H0/(Constants.c*Constants.c*k*k*a*a)*OmegaR0*Theta2;

  //=============================================================================
  // DONE: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Phi,  delta,  v,  ...)
  dPhidx = Psi-ck_Hp*ck_Hp*Phi/3.+ H0*H0/(2.*Hp*Hp)*(OmegaCDM0/a*delta_cdm+Omegab0/a*delta_b+4.*OmegaR0/(a*a)*Theta0);
  ddelta_cdmdx = ck_Hp*v_cdm-3.*dPhidx;
  dv_cdmdx = -v_cdm-ck_Hp*Psi;
  ddelta_bdx = ck_Hp*v_b-3.*dPhidx;
  dv_bdx = -v_b-ck_Hp*Psi+dtau_dx*R*(3.*Theta1+v_b);

  // SET: Photon multipoles (Theta_ell)
  dTheta0_dx = -ck_Hp*Theta1-dPhidx;
  dTheta1_dx = ck_Hp/3.*Theta0-2.*ck_Hp/3.*Theta2+ck_Hp/3.*Psi+ dtau_dx*(Theta1+1./3.*v_b);
  double delta;
  for (int i = 2;i < n_ell_theta-1 ;i++){
    if (i == 2){
      delta = 1.;
    }
    else{
      delta = 0.;
    }
    dydx[Constants.ind_start_theta+i] = i*ck_Hp/(2.*i+1.)*y[Constants.ind_start_theta+i-1]
                                      -(i+1.)*ck_Hp/(2.*i+1.)*y[Constants.ind_start_theta+i+1]
                                      +dtau_dx*(y[Constants.ind_start_theta+i]-1./10.*Theta2*delta);
  }

  dTheta7_dx = ck_Hp*Theta6 - Constants.c*(7.+1.)/(Hp*eta_x)*Theta7 + dtau_dx*Theta7;

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x,  const double k) const{
  return delta_cdm_spline(x, k);
}
double Perturbations::get_delta_b(const double x,  const double k) const{
  return delta_b_spline(x, k);
}
double Perturbations::get_v_cdm(const double x,  const double k) const{
  return v_cdm_spline(x, k);
}
double Perturbations::get_v_b(const double x,  const double k) const{
  return v_b_spline(x, k);
}
double Perturbations::get_Phi(const double x,  const double k) const{
  return Phi_spline(x, k);
}
double Perturbations::get_Psi(const double x,  const double k) const{
  return Psi_spline(x, k);
}
double Perturbations::get_Pi(const double x,  const double k) const{
  return Pi_spline(x, k);
}
double Perturbations::get_Source_T(const double x,  const double k) const{
  return ST_spline(x, k);
}
double Perturbations::get_Source_E(const double x,  const double k) const{
  return SE_spline(x, k);
}
double Perturbations::get_Theta(const double x,  const double k,  const int ell) const{
  return Theta_spline[ell](x, k);
}


//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k,  const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start,  x_end,  npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x, k, 0)   << " ";
    fp << get_Theta(x, k, 1)   << " ";
    fp << get_Theta(x, k, 2)   << " ";
    fp << get_Phi(x, k)       << " ";
    fp << get_Psi(x, k)       << " ";
    fp << get_delta_cdm(x, k) << " ";
    fp << get_delta_b(x, k)   << " ";
    fp << get_v_cdm(x, k)     << " ";
    fp << get_v_b(x, k)       << " ";
    fp << cosmo->eta_of_x(x) << " ";

    //fp << get_Pi(x, k)        << " ";
    //fp << get_Source_T(x, k)  << " ";
    //fp << get_Source_T(x, k) * Utils::j_ell(5,    arg)           << " ";
    //fp << get_Source_T(x, k) * Utils::j_ell(50,   arg)           << " ";
    //fp << get_Source_T(x, k) * Utils::j_ell(500,  arg)           << " ";

    fp << "\n";
  };
  std::for_each(x_array.begin(),  x_array.end(),  print_data);
}
