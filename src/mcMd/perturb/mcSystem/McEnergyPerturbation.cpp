#ifdef  MCMD_PERTURB
#ifndef MCMD_MC_ENERGY_PERTURBATION_CPP
#define MCMD_MC_ENERGY_PERTURBATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McEnergyPerturbation.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <util/ensembles/EnergyEnsemble.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McEnergyPerturbation::McEnergyPerturbation(McSystem& system, 
                                              int size, int rank)
    : LinearPerturbation<McSystem>(system, size, rank)
   {}

   /*
   * Read beta = 1/kT from file
   */
   void McEnergyPerturbation::readParameters(std::istream& in)
   {  
      Perturbation::readParameters(in);
      nParameters_ = Perturbation::getNParameters(); 
   }

   /*
   * Set inverse temperature.
   *
   * \param parameter inverse temperature beta.
   */
   void McEnergyPerturbation::setParameter()
   {  system().energyEnsemble().setTemperature(1.0/parameter_[0]); }

   /*
   * Get inverse temperature of the parent system.
   */
   double McEnergyPerturbation::parameter(int i) const
   {  
      if (i > nParameters_) {
         UTIL_THROW("perturbation parameter index is out of bounds");
      }
      double T = system().energyEnsemble().temperature(); 
      return (1.0/T);
   }

   /*
   * Return the system energy.
   */
   double McEnergyPerturbation::derivative(int i) const
   {  
      if (i > nParameters_) {
         UTIL_THROW("perturbation parameter index is out of bounds");
      }
      return system().potentialEnergy(); 
   }

}

#endif  // ifndef ENERGY_PERTURBATION_CPP
#endif  // ifdef MCMD_PERTURB
