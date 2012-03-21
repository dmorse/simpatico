#ifndef MC_BOND_ENERGY_AVERAGE_CPP
#define MCMD_MC_BOND_ENERGY_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McBondEnergyAverage.h"                    // class header
#include <mcMd/potentials/bond/BondPotential.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   McBondEnergyAverage::McBondEnergyAverage(McSystem& system)
    : AverageDiagnostic<McSystem>(system)
   {}

   /*
   * Evaluate total bond energy.
   */
   void McBondEnergyAverage::sample(long iStep) 
   {  accumulator_.sample(system().bondPotential().energy(), outputFile_); }

}
#endif 
