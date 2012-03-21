#ifndef DDMD_BOND_POTENTIAL_CPP
#define DDMD_BOND_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BondPotential.h"
#include <ddMd/system/System.h>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   BondPotential::BondPotential(System& system)
    : boundaryPtr_(&system.boundary()),
      storagePtr_(&system.bondStorage())
   {}

   /*
   * Constructor (for unit testing).
   */
   BondPotential::BondPotential(Boundary& boundary,
                                GroupStorage<2>& storage)
    : boundaryPtr_(&boundary),
      storagePtr_(&storage)
   {}

   /*
   * Destructor.
   */
   BondPotential::~BondPotential()
   {}

}
#endif
