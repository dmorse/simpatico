/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McPressureAnalyzer.h"

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   McPressureAnalyzer::McPressureAnalyzer(McSystem& system)
    : AverageAnalyzer<McSystem>(system)
   {  setClassName("McPressureAnalyzer"); }

   /* 
   * Evaluate pressure, and add to accumulator.
   */
   void McPressureAnalyzer::compute() 
   {  system().computeStress(value_); }

}
