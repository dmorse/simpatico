#ifndef SIMP_BEND_FORCE_H
#define SIMP_BEND_FORCE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Bend.h"   // base class

#include <cmath>

namespace Simp
{

   using namespace Util;

   /**
   * A BendForce computes derivatives of the angle between two vectors. 
   *
   * Models the angle between a pair of vectors b1 and b2 that are 
   * passed to the functions computeAngle() or computeDerivatives().
   *
   * The scalar cosTheta is the cosine of the angle between b1 and b2.
   *
   * The vectors d1 and d2 are vectors whose elements are the derivatives 
   * of cosTheta with respect to the elements of b1 and b2, respectively.
   * 
   * \ingroup Simp_Interaction_Angle_Module
   */
   struct BendForce : public Bend
   {
   
   public:

      /**
      * Vector of derivatives d1[i] = d(cosTheta)/d(b1[i])
      */
      Vector d1;

      /**
      * Vector of derivatives d2[i] = d(cosTheta)/d(b2[i])
      */
      Vector d2;
   
      /**
      * Compute cosTheta and its derivatives d1 and d2.
      *
      * \param b1  bond vector from atom 0 to 1.
      * \param b2  bond vector from atom 1 to 2.
      */
      void computeDerivatives(const Vector& b1, const Vector& b2);

   };

   // Inline method definitions
   
   /* 
   * Calculate cosTheta and unit vectors.
   */
   inline
   void BendForce::computeDerivatives(const Vector& b1, const Vector& b2)
   {
      Vector u1, u2;

      double b1Abs = b1.abs();
      u1.divide(b1, b1Abs);

      double b2Abs = b2.abs();
      u2.divide(b2, b2Abs);

      cosTheta = u1.dot(u2);

      Vector t;
      t.multiply(u1, cosTheta);
      d1.subtract(u2, t);
      d1 /= b1Abs;

      t.multiply(u2, cosTheta);
      d2.subtract(u1, t);
      d2 /= b2Abs;
   }

} 
#endif
