#ifdef  MCMD_PERTURB
#ifndef HOOMD_MC_PERTURBATION_FACTORY_CPP
#define HOOMD_MC_PERTURBATION_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdMcPerturbationFactory.h"  

// Subclasses of Perturbation
#ifndef SIMP_NOPAIR
#include <mcMd/perturb/mcSystem/McPairPerturbation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <modules/hoomd/potentials/pair/HoomdLJPair.h>
#include <modules/hoomd/potentials/pair/HoomdDpdPair.h>
#endif

#ifdef SIMP_EXTERNAL
#include <mcMd/perturb/mcSystem/McExternalPerturbation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/external/ExternalPotential.h>
#include <modules/hoomd/potentials/external/HoomdPeriodicExternal.h>
#endif

#ifndef SIMP_NOPAIR
#ifdef SIMP_EXTERNAL
#include <mcMd/perturb/mcSystem/McPairExternalPerturbation.h>
#endif
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

      int size = systemPtr_->simulation().communicator().Get_size();
      int rank = systemPtr_->simulation().communicator().Get_rank();

      #ifndef SIMP_NOPAIR
      if (className == "McPairPerturbation") {
         const std::string& interactionClassName = systemPtr_->pairPotential().
            interactionClassName();
         if (interactionClassName == "HoomdLJPair") {
            ptr = new McPairPerturbation<HoomdLJPair> (*systemPtr_, size, rank);
         } else if (interactionClassName == "HoomdDpdPair") {
            ptr = new McPairPerturbation<HoomdDpdPair> (*systemPtr_, size, rank);
         } else { 
            UTIL_THROW("Unsupported pair potential.");
         }
      } 
      #endif
      #ifdef SIMP_EXTERNAL
      else if (className == "McExternalPerturbation") {
         const std::string& interactionClassName = systemPtr_->externalPotential().
            interactionClassName();
         if (interactionClassName == "HoomdPeriodicExternal") {
            ptr = new McExternalPerturbation<HoomdPeriodicExternal> (*systemPtr_, size, rank);
         } else { 
            UTIL_THROW("Unsupported external potential.");
         }
      }
      #endif
      #ifdef SIMP_EXTERNAL
      #ifndef SIMP_NOPAIR
      else if (className == "McPairExternalPerturbation") {
         const std::string& pairInteractionClassName = systemPtr_->pairPotential().
            interactionClassName();
         const std::string& externalInteractionClassName = systemPtr_->externalPotential().
            interactionClassName();
         if (pairInteractionClassName == "HoomdLJPair") {
            if (externalInteractionClassName == "HoomdPeriodicExternal") {
               ptr = new McPairExternalPerturbation<HoomdLJPair,HoomdPeriodicExternal> (*systemPtr_, size, rank);
            } else {
               UTIL_THROW("Unsupported external potential.");
            }
         } else if (pairInteractionClassName == "HoomdDpdPair") {
            if (externalInteractionClassName == "HoomdPeriodicExternal") {
               ptr = new McPairExternalPerturbation<HoomdDpdPair,HoomdPeriodicExternal> (*systemPtr_, size, rank);
            } else {
               UTIL_THROW("Unsupported external potential.");
            }
         } else { 
            UTIL_THROW("Unsupported pair potential.");
         }
      }
      #endif
      #endif

      return ptr;
   }

}

#endif
#endif  // #ifdef  MCMD_PERTURB
