#ifndef MCMD_MC_BOND_ENERGY_AVERAGE_H
#define MCMD_MC_BOND_ENERGY_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/util/AverageDiagnostic.h>  // base class template
#include <mcMd/mcSimulation/McSystem.h>               // base template parameter

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /**
   * McBondEnergyAverage averages of bond potential energy.
   *
   * \ingroup McMd_Diagnostic_Module
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
  
      using AverageDiagnostic<McSystem>::readParameters;
      using AverageDiagnostic<McSystem>::loadParameters;
      using AverageDiagnostic<McSystem>::save;
      using AverageDiagnostic<McSystem>::setup;
      using AverageDiagnostic<McSystem>::output;

   protected:
 
      using AverageDiagnostic<McSystem>::outputFile_;
      using AverageDiagnostic<McSystem>::nSamplePerBlock_;
      using AverageDiagnostic<McSystem>::accumulator_;

   };

}
#endif 
