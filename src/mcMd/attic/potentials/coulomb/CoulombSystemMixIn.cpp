#ifndef MCMD_COULOMB_SYSTEM_MIXIN_CPP
#define MCMD_COULOMB_SYSTEM_MIXIN_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CoulombSystemMixIn.h" 
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>

namespace McMd
{

   /*
   * Constructor.
   */
   CoulombSystemMixIn::CoulombSystemMixIn(System& system)
    : simulationPtr_(&system.simulation()),
      systemPtr_(&system),
      boundaryPtr_(&system.boundary()),
      atomTypesPtr_(&system.simulation().atomTypes())
   {}

} 
#endif
