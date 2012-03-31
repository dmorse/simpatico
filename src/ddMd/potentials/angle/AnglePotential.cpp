#ifndef DDMD_ANGLE_POTENTIAL_CPP
#define DDMD_ANGLE_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AnglePotential.h"
#include <ddMd/simulation/Simulation.h>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   AnglePotential::AnglePotential(Simulation& simulation)
    : boundaryPtr_(&simulation.boundary()),
      storagePtr_(&simulation.angleStorage())
   {}

   /*
   * Constructor (for unit testing).
   */
   AnglePotential::AnglePotential(Boundary& boundary,
                                GroupStorage<3>& storage)
    : boundaryPtr_(&boundary),
      storagePtr_(&storage)
   {}

   /*
   * Destructor.
   */
   AnglePotential::~AnglePotential()
   {}

}
#endif
