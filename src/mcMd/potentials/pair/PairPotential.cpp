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
      if (!energy_.isSet()) {
         computeEnergy();
      }
      return energy_.value();
   }

   /*
   * Mark the stress as unknown.
   */
   void PairPotential::unsetStress()
   {  stress_.unset(); }

   /*
   * Get the nonbonded stress tensor.
   */
   void PairPotential::computeStress(Tensor& stress)
   {
      // If necessary, compute stress 
      if (!stress_.isSet()) {
         computeStress();
      }

      // Get pair stress tensor
      stress = stress_.value();
   }

   /*
   * Get the nonbonded x, y, z pressures
   */
   void PairPotential::computeStress(Vector& pressures)
   {
      // If necessary, compute stress 
      if (!stress_.isSet()) {
         computeStress();
      }

      // Get diagonal components of pair stress.
      for (int i=0; i < Dimension; ++i) {
         pressures[i] = stress_.value()(i, i);
      }
   }

   /*
   * Get the nonbonded pressure
   */
   void PairPotential::computeStress(double& pressure)
   {
      // If necessary, compute stress 
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
