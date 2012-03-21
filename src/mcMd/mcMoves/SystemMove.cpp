#ifndef SYSTEM_MOVE_CPP
#define MCMD_SYSTEM_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   SystemMove::SystemMove(McSystem& system) :
      McMove(system.simulation()),
      systemPtr_(&system),
      boundaryPtr_(&system.boundary())
   {
      isothermalPtr_ = &system.energyEnsemble();
   }

   /* 
   * Destructor.
   */
   SystemMove::~SystemMove() 
   {}

}
#endif
