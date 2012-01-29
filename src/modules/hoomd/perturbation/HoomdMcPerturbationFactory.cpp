#ifdef  MCMD_PERTURB
#ifndef HOOMD_MC_PERTURBATION_FACTORY_CPP
#define HOOMD_MC_PERTURBATION_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdMcPerturbationFactory.h"  

// Subclasses of Perturbation
#ifndef MCMD_NOPAIR
#include <mcMd/perturb/mcSystem/McPairPerturbation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotential.h>

#include <modules/hoomd/potentials/pair/HoomdLJPair.h>
#include <modules/hoomd/potentials/pair/HoomdDpdPair.h>
#endif

namespace McMd
{

   using namespace Util;

   HoomdMcPerturbationFactory::HoomdMcPerturbationFactory(McSystem& system)
    : systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Peturbation subclass className.
   */
   Perturbation* HoomdMcPerturbationFactory::factory(const std::string &className) const
   {
      Perturbation *ptr = 0;

      #ifndef MCMD_NOPAIR
      if (className == "McPairPerturbation") {
         const std::string& evaluatorClassName = systemPtr_->pairPotential().
            evaluatorClassName();
         if (evaluatorClassName == "HoomdLJPair") {
            ptr = new McPairPerturbation<HoomdLJPair> (*systemPtr_);
         } else if (evaluatorClassName == "HoomdDpdPair") {
            ptr = new McPairPerturbation<HoomdDpdPair> (*systemPtr_);
         } 
      } 
      #endif
      return ptr;
   }

}

#endif
#endif  // #ifdef  MCMD_PERTURB
