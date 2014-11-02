/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
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
