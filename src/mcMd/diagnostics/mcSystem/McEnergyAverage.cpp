#ifndef MCMD_MC_ENERGY_AVERAGE_CPP
#define MCMD_MC_ENERGY_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
    : AverageDiagnostic<McSystem>(system)
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
#endif 
