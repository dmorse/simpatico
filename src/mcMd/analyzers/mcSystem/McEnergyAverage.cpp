/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McEnergyAverage.h"                        // class header
//#include <util/misc/FileMaster.h>  
//#include <util/archives/Serializable_includes.h>


#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   McEnergyAverage::McEnergyAverage(McSystem& system)
    : AverageAnalyzer<McSystem>(system)
   {  setClassName("McEnergyAverage"); }


   /* 
   * Evaluate energy, and add to accumulator.
   */
   void McEnergyAverage::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         accumulator_.sample(system().potentialEnergy(), outputFile_);
      }
   }
   
}
