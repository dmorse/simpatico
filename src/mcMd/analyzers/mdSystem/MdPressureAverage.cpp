#ifndef MCMD_MD_PRESSURE_AVERAGE_CPP
#define MCMD_MD_PRESSURE_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdPressureAverage.h"        // class header

namespace McMd
{

   /* 
   * Constructor.
   */
   MdPressureAverage::MdPressureAverage(MdSystem& system)
    : PressureAverage<MdSystem>(system)
   {  setClassName("MdPressureAverage"); }

}
#endif
