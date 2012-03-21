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

   class System;

   /*
   * Constructor.
   */
   Integrator::Integrator(System& system)
     : systemPtr_(&system)
   {}

   /*
   * Destructor.
   */
   Integrator::~Integrator()
   {}

}
#endif
