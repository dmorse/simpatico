#ifndef MCMD_MC_PRESSURE_AVERAGE_CPP
#define MCMD_MC_PRESSURE_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/


#include "McPressureAverage.h"        // class header

namespace McMd
{

   /* 
   * Constructor.
   */
   McPressureAverage::McPressureAverage(McSystem& system)
    : PressureAverage<McSystem>(system)
   {}

}
#endif
