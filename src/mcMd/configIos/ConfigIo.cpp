#ifndef CONFIG_IO_CPP
#define CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigIo.h"
#include <mcMd/simulation/System.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor. 
   */
   ConfigIo::ConfigIo(System& system)
    : boundaryPtr_(&system.boundary()),
      systemPtr_(&system),
      simulationPtr_(&system.simulation())
   {}

   /* 
   * Destructor.   
   */
   ConfigIo::~ConfigIo() 
   {}

} 
#endif
