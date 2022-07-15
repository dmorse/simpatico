/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
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
