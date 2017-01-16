#include "PairPotential.h"
#include <util/space/Vector.h>

namespace McMd 
{

   using namespace Util;   

   /*
   * Mark the energy as unknown.
   */
   void PairPotential::unsetEnergy()
   {  energy_.unset(); }

   /**
   * Return the precalculated total pair energy for this System.
   */
   double PairPotential::energy()
   {
      computeEnergy();
      return energy_.value();
   }

   /*
   * Mark the stress as unknown.
   */
   void PairPotential::unsetStress()
   {  stress_.unset(); }

   /*
   * Compute x, y, z nonbonded pressures.
   */
   void PairPotential::computeStress(Tensor& stress)
   {
      computeStress();
      stress = stress_.value();
   }

   /*
   * Compute the total nonbonded pressure
   */
   void PairPotential::computeStress(Vector& pressures)
   {
      // Compute stress if necessary.
      computeStress();

      // Get diagonal components of pair stress.
      for (int i=0; i < Dimension; ++i) {
         pressures[i] = stress_.value()(i, i);
      }
   }

   /*
   * Compute the total nonbonded pressure
   */
   void PairPotential::computeStress(double& pressure)
   {
      // Compute stress if necessary.
      computeStress();

      // Compute pressure = average of diagonal components.
      pressure = 0.0;
      for (int i=0; i < Dimension; ++i) {
         pressure += stress_.value()(i, i);
      }
      pressure = pressure/double(Dimension);
   }

}
