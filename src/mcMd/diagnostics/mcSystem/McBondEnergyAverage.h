#ifndef MCMD_MC_BOND_ENERGY_AVERAGE_H
#define MCMD_MC_BOND_ENERGY_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/util/AverageDiagnostic.h>  // base class template
#include <mcMd/mcSimulation/McSystem.h>               // base template parameter

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /**
   * McBondEnergyAverage averages of total potential energy.
   *
   * \ingroup Diagnostic_Module
   */
   class McBondEnergyAverage : public AverageDiagnostic<McSystem>
   {
   
   public:

      /**   
      * Constructor.
      */
      McBondEnergyAverage(McSystem& system);

      /**
      * Evaluate total bond energy, and add to ensemble. 
      */
      virtual void sample(long iStep);
  
      using AverageDiagnostic<McSystem>::readParam;
      using AverageDiagnostic<McSystem>::output;

   protected:
 
      using AverageDiagnostic<McSystem>::outputFile_;
      using AverageDiagnostic<McSystem>::accumulator_;
      using AverageDiagnostic<McSystem>::nSamplePerBlock_;

   };

}
#endif 
