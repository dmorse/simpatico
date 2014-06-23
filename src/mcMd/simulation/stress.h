#ifndef MCMD_STRESS_H
#define MCMD_STRESS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>               // member template
#include <util/space/Tensor.h>               // member template

namespace McMd
{

   using namespace Util;

   /*
   * The functions defined in this file are used by function templates that calculate
   * the isotropic pressure (double), x-y-z diagonal stress components (Vector),
   * and the full stress tensor (Tensor) using a common template. 
   *
   * The overloaded incrementPairStress correctly takes care of the addition of the 
   * contribution from the force between a pair of particles to either a double
   * precision pressure or to a Vector or Tensor of stress components. 
   *
   * The overloaded normalizeStress function divides the trace of the stress tensor 
   * by 3.0 (the dimensionality of space) for a double precision pressure, but does
   * nothing when called with a Vector or Tensor of stress components.
   */

   /*
   * Add a a pair contribution to isotropic pressure.
   */
   inline void incrementPairStress(const Vector& f, const Vector& dr, double& pressure)
   {  pressure += f.dot(dr); }

   /*
   * Add a a pair contribution to x, y, and z stress components.
   */
   inline void incrementPairStress(const Vector& f, const Vector& dr, Vector& stress)
   {
      for (int i = 0; i < Dimension; ++i) {
         stress[i] += f[i]*dr[i];
      }
   }

   /*
   * Add a pair contribution to the stress tensor.
   */
   inline void incrementPairStress(const Vector& f, const Vector& dr, Tensor& stress)
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            stress(i, j) += f[i]*dr[j];
         }
      }
   }

   /*
   * Divide trace of stress by three to obtain pressure.
   */
   inline void normalizeStress(double& pressure)
   {  pressure /= 3.0; }

   /*
   * Do nothing
   */
   inline void normalizeStress(Vector& stress)
   { }

   /*
   * Do nothing
   */
   inline void normalizeStress(Tensor& stress)
   { }

}
#endif
