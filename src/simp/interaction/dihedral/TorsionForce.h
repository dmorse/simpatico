#ifndef SIMP_TORSION_FORCE_H
#define SIMP_TORSION_FORCE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Torsion.h"

#include <cmath>

namespace Simp
{

   using namespace Util;

   /**
   * Computes derivatives of dihedral angle with respect to bond vectors.
   *
   * Models the dihedral angle formed by three sequentil bond vectors b1, b2,
   * and b3. See \ref Simp_Interaction_Dihedral_Module in the file dihedral.mod for the 
   * definition of the dihedral angle phi and its relationship to these 
   * bond vectors
   *
   * These 3 bond vectors must be passed to the method computeDerivatives()
   * which computes the cosine of the dihedral angle and the derivatives of
   * the cosine with respect to the elements of the bond vectors. This 
   * function stores its results in the public members cosPhi, d1, d2, and
   * d3. Upon return:
   *
   * The scalar member cosPhi is the cosine of the dihdral angle phi.
   *
   * The elements of the vectors d1, d2, and d3 are derivatives of cosPhi 
   * with respect to the elements of the bond vectors b1, b2, and b3.
   *
   * The function computeAngle(), which is inherited from the Torsion
   * base class, computes only cosPhi, but not the derivatives. 
   *
   * \ingroup Simp_Interaction_Dihedral_Module
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
      *
      * \return 0 if normal, 1 for divide by zero error
      */
      bool computeDerivatives(const Vector& b1, const Vector& b2, const Vector& b3);
 
   };

   // Inline method definitions

   /* 
   * Calculate cosPhi.
   */ 
   inline bool 
   TorsionForce::computeDerivatives(
                    const Vector& b1, const Vector& b2, const Vector& b3)
   {
      Vector u1, u2, t1, t2;
      double r1, r2;

      u1.cross(b1, b2);
      r1 = u1.square();
      if (r1 < 1.0E-10) {
         return 1; // Error code
      }
      r1 = sqrt(r1);
      u1 /= r1;

      u2.cross(b2, b3);
      r2 = u2.square();
      if (r2 < 1.0E-10) {
         return 1; // Error code
      }
      r2 = sqrt(r2);
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

      return 0; // Normal completion
   }

} 
#endif
