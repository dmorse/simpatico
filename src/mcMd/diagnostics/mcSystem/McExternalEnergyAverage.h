#ifdef INTER_EXTERNAL
#ifndef MCMD_MC_EXTERNAL_ENERGY_AVERAGE_H
#define MCMD_MC_EXTERNAL_ENERGY_AVERAGE_H

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
   * McExternalEnergyAverage averages of total external energy.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class McExternalEnergyAverage : public AverageDiagnostic<McSystem>
   {
   
   public:

      /**   
      * Constructor.
      *
      * \param system parent McSystem
      */
      McExternalEnergyAverage(McSystem& system);

      /* 
      * Evaluate external energy per particle, and add to ensemble. 
      *
      * \param iStep step counter index
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
#endif 
