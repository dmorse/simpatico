#ifdef  MCMD_PERTURB
#ifndef MC_PERTURBATION_FACTORY_CPP
#define MC_PERTURBATION_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McPerturbationFactory.h"  

// Subclasses of Perturbation
#include "McEnergyPerturbation.h"
#ifndef MCMD_NOPAIR
#include "McPairPerturbation.h"
#include <mcMd/potentials/pair/LJPair.h>
#include <mcMd/potentials/pair/DpdPair.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#ifdef MCMD_EXTERNAL
#include "McExternalPerturbation.h"
#endif
#ifndef MCMD_NOPAIR
#ifdef MCMD_EXTERNAL
#include "McPairExternalPerturbation.h"
#endif
#endif

namespace McMd
{

   using namespace Util;

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
      #ifdef MCMD_EXTERNAL 
      if (className == "McExternalPerturbation") {
         ptr = new McExternalPerturbation(*systemPtr_);
      } 
      #endif
      #ifndef MCMD_NOPAIR
      else if (className == "McPairPerturbation") {
         const std::string& evaluatorClassName = systemPtr_->pairPotential().
            evaluatorClassName();
         if (evaluatorClassName == "LJPair") {
            ptr = new McPairPerturbation<LJPair> (*systemPtr_);
         } else if (evaluatorClassName == "DpdPair") {
            ptr = new McPairPerturbation<DpdPair> (*systemPtr_);
         } else {
            UTIL_THROW("Unsupported pair potential.");
         }
      } 
      #endif
      #ifndef MCMD_NOPAIR
      #ifdef MCMD_EXTERNAL
      else if (className == "McPairExternalPerturbation") {
         const std::string& evaluatorClassName = systemPtr_->pairPotential().
            evaluatorClassName();
         if (evaluatorClassName == "LJPair") {
            ptr = new McPairExternalPerturbation<LJPair> (*systemPtr_);
         } else if (evaluatorClassName == "DpdPair") {
            ptr = new McPairExternalPerturbation<DpdPair> (*systemPtr_);
         } else {
            UTIL_THROW("Unsupported pair potential.");
         }
      }
      #endif 
      #endif
      /*#ifdef MCMD_EXTERNAL
      else if (className == "McExternalPerturbation") {
         const std::string& evaluatorClassName = systemPtr_->externalPotential().
            evaluatorClassName();
         if (evaluatorClassName == "TanhCosineExternal") {
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
