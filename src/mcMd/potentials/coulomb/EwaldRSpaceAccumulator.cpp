
#include "EwaldRSpaceAccumulator.h"
#include <mcMd/potentials/pair/PairPotential.h>

namespace McMd
{

   using namespace Util;

   /*
   * Return the r-space energy.
   */
   double EwaldRSpaceAccumulator::rSpaceEnergy() 
   {
      if (!rSpaceEnergy_.isSet()) {
         UTIL_CHECK(pairPotentialPtr_);
         pairPotentialPtr_->computeEnergy();
      } 
      return rSpaceEnergy_.value(); 
   }

   /*
   * Return the r-space stress.
   */
   Tensor EwaldRSpaceAccumulator::rSpaceStress()
   { 
      if (!rSpaceStress_.isSet()) {
         UTIL_CHECK(pairPotentialPtr_);
         pairPotentialPtr_->computeStress();
      }
      return rSpaceStress_.value(); 
   }
     
}
