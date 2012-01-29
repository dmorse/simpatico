#ifdef  MCMD_EXTERNAL
#ifndef EXTERNAL_FACTORY_CPP
#define EXTERNAL_FACTORY_CPP

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

#include <mcMd/potentials/external/BoxExternal.h>
#include <mcMd/potentials/external/OrthoBoxExternal.h>
#include <mcMd/potentials/external/SlitExternal.h>
#include <mcMd/potentials/external/TanhCosineExternal.h>

namespace McMd
{

   /**
   * Default constructor.
   */
   ExternalFactory::ExternalFactory(System& system)
    : systemPtr_(&system)
   {}

   /*
   * Return a pointer to a new ExternalPotential, if possible.
   */
   ExternalPotential* ExternalFactory::factory(const std::string& name) const
   {
      ExternalPotential* ptr = 0;

      ptr = trySubfactories(name);
      if (ptr) return ptr;

      if (name == "BoxExternal") {
         ptr = new ExternalPotentialImpl<BoxExternal>(*systemPtr_);
      } else
      if (name == "OrthoBoxExternal") {
         ptr = new ExternalPotentialImpl<OrthoBoxExternal>(*systemPtr_);
      } else
      if (name == "SlitExternal") {
         ptr = new ExternalPotentialImpl<SlitExternal>(*systemPtr_);
      } else
      if (name == "TanhCosineExternal") {
         ptr = new ExternalPotentialImpl<TanhCosineExternal>(*systemPtr_);
      }
      return ptr;
   }

}
#endif
#endif
