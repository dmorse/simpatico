#ifdef MCMD_EXTERNAL
#ifndef HOOMD_EXTERNAL_FACTORY_CPP
#define HOOMD_EXTERNAL_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/external/ExternalFactory.h>

#include <mcMd/simulation/System.h>

// ExternalPotential interfaces and implementation classes
#include <mcMd/potentials/external/ExternalPotential.h>
#include <mcMd/potentials/external/ExternalPotentialImpl.h>

// External Potential evaluator classes
#include <inter/external/TanhCosineExternal.h>

#include <modules/hoomd/potentials/external/HoomdLamellarExternal.h>


#include <hoomd/HOOMDMath.h>
#include <hoomd/EvaluatorExternalLamellar.h>
#include <hoomd/AllDriverPotentialExternalGPU.cuh>

#include "HoomdExternalFactory.h"

namespace McMd
{

   /**
   * Default constructor.
   */
   HoomdExternalFactory::HoomdExternalFactory(System &system)
      : ExternalFactory(system), systemPtr_(&system) 
   {}

   /*
   * Return a pointer to a new ExternalInteration, if possible.
   */
   ExternalPotential* 
   HoomdExternalFactory::factory(const std::string& name) const
   {
      ExternalPotential* ptr = 0;
      if (name == "HoomdLamellarExternal") {
         ptr = new ExternalPotentialImpl< HoomdLamellarExternal> (*systemPtr_);
      }
      return ptr;
   }

   /**
   * return a ForceCompute
   */
   boost::shared_ptr<ForceCompute> HoomdExternalFactory::hoomdFactory(
      ExternalPotential& potential,
      System &system,
      boost::shared_ptr<SystemDefinition> systemDefinitionSPtr)
   {
      std::string className = potential.evaluatorClassName();
      boost::shared_ptr<ForceCompute> externalSPtr;

      if (className == "HoomdLamellarExternal") {
         externalSPtr = hoomdFactoryImpl<EvaluatorExternalLamellar, gpu_compute_lamellar_forces, HoomdLamellarExternal >(&potential, system, systemDefinitionSPtr );
      } else
         UTIL_THROW("Unsupported Hoomd potential." );

      return externalSPtr; 
   }


}
#endif
#endif
