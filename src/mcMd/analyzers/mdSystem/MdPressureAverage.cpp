/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
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
