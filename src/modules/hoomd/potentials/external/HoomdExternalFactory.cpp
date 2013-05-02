#ifdef INTER_EXTERNAL
#ifndef EXTERNAL_FACTORY_CPP
#define EXTERNAL_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/external/ExternalFactory.h>

#include <mcMd/simulation/System.h>

// ExternalPotential interfaces and implementation classes
#include <mcMd/potentials/external/ExternalPotential.h>
#include <mcMd/potentials/external/ExternalPotentialImpl.h>

// External Potential evaluator classes
#include <inter/external/LamellarOrderingExternal.h>
#include <inter/external/LocalLamellarOrderingExternal.h>

#include <modules/hoomd/potentials/external/HoomdLamellarOrderingExternal.h>
#include <modules/hoomd/potentials/external/HoomdLocalLamellarOrderingExternal.h>

#include <hoomd/HOOMDMath.h>
#include <hoomd/LocalExternalParams.h>
#include <hoomd/EvaluatorExternalPeriodic.h>
#include <hoomd/EvaluatorLocalExternal.h>
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
      if (name == classNameHoomdLamellarOrdering) {
         ptr = new ExternalPotentialImpl<HoomdLamellarOrderingExternal> (*systemPtr_);
      } else
      if (name == classNameHoomdLocalLamellarOrdering) {
         ptr = new ExternalPotentialImpl<HoomdLocalLamellarOrderingExternal> (*systemPtr_);
      }
      return ptr;
   }

   /** 
   * return a HoomdExternalPotential
   */
   HoomdExternalPotential *HoomdExternalFactory::hoomdExternalPotentialPtr(ExternalPotential&
      potential)
   {
      std::string name = potential.interactionClassName();
      HoomdExternalPotential* ptr = 0;

      if (name == classNameHoomdLamellarOrdering) {
         ptr = dynamic_cast< HoomdExternalPotential *>(
            dynamic_cast< HoomdLamellarOrderingExternal * >(
            &(dynamic_cast< ExternalPotentialImpl< HoomdLamellarOrderingExternal > * >
            (&potential))->interaction()));
      } else
      if (name == classNameHoomdLocalLamellarOrdering) {
         ptr = dynamic_cast< HoomdExternalPotential *>(
            dynamic_cast< HoomdLocalLamellarOrderingExternal * >(
            &(dynamic_cast< ExternalPotentialImpl< HoomdLocalLamellarOrderingExternal > * >
            (&potential))->interaction()));
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
      // first get a pointer to a HoomdExternalPotential
      HoomdExternalPotential* ptr = hoomdExternalPotentialPtr(potential);
             
      std::string className = potential.interactionClassName();
      boost::shared_ptr<ForceCompute> externalSPtr;

      if (className == classNameHoomdLamellarOrdering) {
         externalSPtr = hoomdFactoryImpl<EvaluatorExternalPeriodic, gpu_compute_periodic_forces, 
                                         classNameHoomdLamellarOrdering >(ptr, system, systemDefinitionSPtr );
      } else 
      if (className == classNameHoomdLocalLamellarOrdering) {
         externalSPtr = hoomdFactoryImpl<EvaluatorLocalExternal, gpu_compute_local_forces, 
                                         classNameHoomdLocalLamellarOrdering >(ptr, system, systemDefinitionSPtr );
      } else 
         UTIL_THROW("Unsupported Hoomd potential." );
      
      return externalSPtr; 
   }


}
#endif
#endif
