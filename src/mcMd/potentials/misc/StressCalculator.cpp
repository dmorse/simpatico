/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StressCalculator.h"
#include <util/space/Vector.h>

namespace McMd 
{

   using namespace Util;   

   /*
   * Mark the stress as unknown.
   */
   void StressCalculator::unsetStress()
   {  stress_.unset(); }

   /*
   * Get the nonbonded stress tensor.
   */
   void StressCalculator::computeStress(Tensor& stress)
   {
      // If necessary, compute stress tensor
      if (!stress_.isSet()) {
         computeStress();
      }

      // Get full stress tensor
      stress = stress_.value();
   }

   /*
   * Get the nonbonded x, y, z pressures
   */
   void StressCalculator::computeStress(Vector& pressures)
   {
      // If necessary, compute stress tensor
      if (!stress_.isSet()) {
         computeStress();
      }

      // Get diagonal components of stress tensor (pressures)
      for (int i=0; i < Dimension; ++i) {
         pressures[i] = stress_.value()(i, i);
      }
   }

   /*
   * Get the isotropic pressure = Tr(stress)/3
   */
   void StressCalculator::computeStress(double& pressure)
   {
      // If necessary, compute stress tensor
      if (!stress_.isSet()) {
         computeStress();
      }

      // Get pressure = average of diagonal components.
      pressure = 0.0;
      for (int i=0; i < Dimension; ++i) {
         pressure += stress_.value()(i, i);
      }
      pressure = pressure/double(Dimension);
   }

}
