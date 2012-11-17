#ifndef DDMD_EXTERNAL_FACTORY_CPP
#define DDMD_EXTERNAL_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/potentials/external/ExternalFactory.h>
#include <ddMd/simulation/Simulation.h>

// ExternalPotential interface and implementation classes
#include <ddMd/potentials/external/ExternalPotential.h>
#include <ddMd/potentials/external/ExternalPotentialImpl.h>

// External interaction classes
#include <inter/external/TanhCosineExternal.h>

namespace DdMd
{

   using namespace Inter;

   /**
   * Default constructor.
   */
   ExternalFactory::ExternalFactory(Simulation& simulation)
    : Factory<ExternalPotential>(),
      simulationPtr_(&simulation)
   {}

   /*
   * Return a pointer to a new ExternalPotential, if possible.
   */
   ExternalPotential* 
   ExternalFactory::factory(const std::string& name) const
   {
      ExternalPotential* ptr = 0;

      // Try subfactories first.
      ptr = trySubfactories(name);
      if (ptr) return ptr;

      if (name == "TanhCosineExternal") {
         ptr = new ExternalPotentialImpl<TanhCosineExternal>(*simulationPtr_);
      } 
      return ptr;
   }

}
#endif
