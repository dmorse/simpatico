#ifndef DDMD_INTEGRATOR_CPP
#define DDMD_INTEGRATOR_CPP

#include "Integrator.h"

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   class Simulation;

   /*
   * Constructor.
   */
   Integrator::Integrator(Simulation& simulation)
     : simulationPtr_(&simulation)
   {}

   /*
   * Destructor.
   */
   Integrator::~Integrator()
   {}

}
#endif
