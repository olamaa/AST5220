#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // Recombination parameters
  double Yp          = 0.0;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  //cosmo.info();
  
  // Output background evolution quantities
  //cosmo.output("cosmology.txt");

  // Do the supernova fits. Uncomment when you are ready to run this
  // Make sure you read the comments on the top of src/SupernovaFitting.h
  
  
  //mcmc_fit_to_supernova_data("data/supernovadata.txt", "results_supernovafitting.txt"); // PROJECT 1


  // Remove when module is completed
  //return 0;

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  
  // Remove when module is completed
  //return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  //std::cout << "line 73" << std::endl;
  Perturbations pert(&cosmo, &rec);
  //std::cout << "line 75" << std::endl;
  pert.solve();
  //std::cout << "line 77" << std::endl;
  pert.info();
  //std::cout << "line 79" << std::endl;
  
  // Output perturbation quantities
  double kvalue0 = 0.001 / Constants.Mpc;
  double kvalue1 = 0.01 / Constants.Mpc;
  double kvalue2 = 0.1 / Constants.Mpc;
  double kvalue3 = 1. / Constants.Mpc;
  pert.output(kvalue0, "perturbations_k_0_001.txt");
  pert.output(kvalue1, "perturbations_k_0_01.txt");
  pert.output(kvalue2, "perturbations_k_0_1.txt");
  pert.output(kvalue3, "perturbations_k_1.txt");
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
