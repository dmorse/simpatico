#ifndef MCMD_MC_STRESS_AUTOCORRELATION_CPP
#define MCMD_MC_STRESS_AUTOCORRELATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/


#include "McStressAutoCorrelation.h"        // class header

namespace McMd
{

   /* 
   * Constructor.
   */
   McStressAutoCorrelation::McStressAutoCorrelation(McSystem& system)
    : StressAutoCorrelation<McSystem>(system)
   {  setClassName("McStressAutoCorrelation"); }

   /* 
   * Destructor.
   */
   McStressAutoCorrelation::~McStressAutoCorrelation()
   {}

   /* 
   * Evaluate pressure, and add to accumulator.
   */
   void McStressAutoCorrelation::sample(long iStep)
   {
      double pressure;
      double temperature;
      McSystem& sys=system(); 
      sys.computeStress(pressure);

      DArray<double> elements;
      elements.allocate(9);

      Tensor total;

      if (isAtInterval(iStep)){
         sys.computeVirialStress(total);
         temperature = sys.energyEnsemble().temperature();

         elements[0] = (total(0,0) - pressure / 3.0) / (10.0 * temperature);
         elements[1] = (total(0,1) + total(1,0)) / 2.0 / (10.0 * temperature);
         elements[2] = (total(0,2) + total(2,0)) / 2.0 / (10.0 * temperature);
         elements[3] = elements[1];
         elements[4] = (total(1,1) - pressure / 3.0) / (10.0 * temperature);
         elements[5] = (total(1,2) + total(2,1)) / 2.0 / (10.0 * temperature);
         elements[6] = elements[2];
         elements[7] = elements[5];
         elements[8] = (total(2,2) - pressure / 3.0) / (10.0 * temperature);

         accumulator_.sample(elements);
     }
   }

}
#endif
