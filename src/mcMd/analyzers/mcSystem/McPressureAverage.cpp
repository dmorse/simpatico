#ifndef MCMD_MC_PRESSURE_AVERAGE_CPP
#define MCMD_MC_PRESSURE_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   {  setClassName("McPressureAverage"); }

   /* 
   * Destructor.
   */
   McPressureAverage::~McPressureAverage()
   {}

}
#endif
