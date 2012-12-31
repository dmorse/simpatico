#ifndef INTER_TORSION_FORCE_H
#define INTER_TORSION_FORCE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Torsion.h"

#include <cmath>

namespace Inter
{

   using namespace Util;

   /**
   * Calculates derivatives of torsion angle involving 3 bonds.
   *
   * Models the dihedral angle of vectors b1, b2, and b3 that are 
   * passed to the function computeAngle() or computeDerivatives().
   *
   * The scalar cosPhi is the cosine of the dihdral angle phi.
   *
   * The elements of the vectors d1, d2, and d3 (if calculated) are
   * derivatives of cosPhi with respect to the elements of b1, b2, 
   * and b3, respectively. 
   * 
   * \ingroup Inter_Dihedral_Module
   */
   struct TorsionForce : public Torsion
   {
   
   public:

      /**
      * Vector of derivatives d1[i] = d(cosPhi)/d(b1[i])
      */
      Vector d1;

      /**
      * Vector of derivatives d2[i] = d(cosPhi)/d(b2[i])
      */
      Vector d2;
   
      /**
      * Vector of derivatives d3[i] = d(cosPhi)/d(b3[i])
      */
      Vector d3;
   
      /**
      * Compute cosPhi and derivatives.
      *
      * \param b1  bond vector from atom 0 to 1.
      * \param b2  bond vector from atom 1 to 2.
      * \param b3  bond vector from atom 2 to 3.
      */
      void computeDerivatives(const Vector& b1, const Vector& b2, const Vector& b3);
 
   };

   // Inline method definitions

   /* 
   * Calculate cosPhi.
   */ 
   inline void 
   TorsionForce::computeDerivatives(
                    const Vector& b1, const Vector& b2, const Vector& b3)
   {
      Vector u1, u2, t1, t2;
      double r1, r2;

      u1.cross(b1, b2);
      r1 = u1.abs();
      u1 /= r1;

      u2.cross(b2, b3);
      r2 = u2.abs();
      u2 /= r2;

      cosPhi = u1.dot(u2);

      t1.multiply(u1, -cosPhi);
      t1 += u2;
      t1 *= 1.0 / r1;

      t2.multiply(u2, -cosPhi);
      t2 += u1;
      t2 *= 1.0 / r2;

      d1.cross(b2, t1);
      d3.cross(t2, b2);

      d2.cross(t1, b1);
      t1.cross(b3, t2);
      d2 += t1;
   }

} 
#endif
