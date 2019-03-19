/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McEnergyAverage.h"                

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McEnergyAverage::McEnergyAverage(McSystem& system)
    : AverageAnalyzer<McSystem>(system)
   {  setClassName("McEnergyAverage"); }

   void McEnergyAverage::compute() {
      value_ = system().potentialEnergy();
   }

}
