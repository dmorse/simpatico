/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EnergyCalculator.h"

namespace McMd 
{

   using namespace Util;   

   /*
   * Mark the energy as unknown (default implementation).
   */
   void EnergyCalculator::unsetEnergy()
   {  energy_.unset(); }

   /**
   * Return the precalculated total pair energy for this System.
   */
   double EnergyCalculator::energy()
   {
      if (!energy_.isSet()) {
         computeEnergy();
      }
      return energy_.value();
   }

}
