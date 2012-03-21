#ifdef  MCMD_PERTURB
#ifndef MC_PERTURBATION_FACTORY_CPP
#define MCMD_MC_PERTURBATION_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McPerturbationFactory.h"  

// Subclasses of Perturbation
#include "McEnergyPerturbation.h"
#ifndef INTER_NOPAIR
#include "McPairPerturbation.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <inter/pair/LJPair.h>
#include <inter/pair/DpdPair.h>
#endif

#ifdef INTER_EXTERNAL
#include "McExternalPerturbation.h"
#ifndef INTER_NOPAIR
#include "McPairExternalPerturbation.h"
#endif
#endif

namespace McMd
{

   using namespace Util;
   using namespace Inter;

   McPerturbationFactory::McPerturbationFactory(McSystem& system)
    : systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Peturbation subclass className.
   */
   Perturbation* McPerturbationFactory::factory(const std::string &className) const
   {
      Perturbation *ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "McEnergyPerturbation") {
         ptr = new McEnergyPerturbation(*systemPtr_);
      }
      #ifdef INTER_EXTERNAL 
      if (className == "McExternalPerturbation") {
         ptr = new McExternalPerturbation(*systemPtr_);
      } 
      #endif
      #ifndef INTER_NOPAIR
      else if (className == "McPairPerturbation") {
         const std::string& interactionClassName = systemPtr_->pairPotential().
            interactionClassName();
         if (interactionClassName == "LJPair") {
            ptr = new McPairPerturbation<LJPair> (*systemPtr_);
         } else if (interactionClassName == "DpdPair") {
            ptr = new McPairPerturbation<DpdPair> (*systemPtr_);
         } else {
            UTIL_THROW("Unsupported pair potential.");
         }
      } 
      #endif
      #ifndef INTER_NOPAIR
      #ifdef INTER_EXTERNAL
      else if (className == "McPairExternalPerturbation") {
         const std::string& interactionClassName = systemPtr_->pairPotential().
            interactionClassName();
         if (interactionClassName == "LJPair") {
            ptr = new McPairExternalPerturbation<LJPair> (*systemPtr_);
         } else if (interactionClassName == "DpdPair") {
            ptr = new McPairExternalPerturbation<DpdPair> (*systemPtr_);
         } else {
            UTIL_THROW("Unsupported pair potential.");
         }
      }
      #endif 
      #endif
      /*#ifdef INTER_EXTERNAL
      else if (className == "McExternalPerturbation") {
         const std::string& interactionClassName = systemPtr_->externalPotential().
            interactionClassName();
         if (interactionClassName == "TanhCosineExternal") {
            ptr = new McExternalPerturbation<TanhCosineExternal> (*systemPtr_);
         } else
            UTIL_THROW("Unsupported external potential.")
      } 
      #endif*/
      return ptr;
   }

}

#endif
#endif  // #ifdef  MCMD_PERTURB
