#ifndef SIMP_BEND_H
#define SIMP_BEND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>

#include <cmath>

namespace Simp
{

   using namespace Util;

   /**
   * A Bend calculates the angle between two vectors.
   *
   * Models the angle between a pair of vectors b1 and b2 that are 
   * passed to the function computeAngle().
   *
   * The scalar cosTheta is the cosine of the angle between b1 and b2.
   * 
   * \ingroup Simp_Interaction_Angle_Module
   */
   struct Bend
   {
   
   public:

      /**
      * Cosine of angle between vectors b1 and b2.
      */
      double cosTheta;

      /**
      * Compute cosTheta.
      *
      * \param b1  bond vector from atom 0 to 1.
      * \param b2  bond vector from atom 1 to 2.
      */
      void computeAngle(const Vector& b1, const Vector& b2);
 
      /**
      * Return value of sin(theta) for precomputed cos(theta).
      */
      double sinTheta() const;

      /**
      * Return value of theta in radians for precomputed cos(theta).
      *
      * Returns principal value, in range 0 < theta < Pi.
      */
      double theta() const;

   };

   // Inline method definitions
   
   /* 
   * Calculate unit bond vectors and cosTheta.
   */
   inline
   void Bend::computeAngle(const Vector& b1, const Vector& b2) 
   {
      double d = sqrt(b1.square()*b2.square());
      cosTheta = b1.dot(b2)/d;
      if (cosTheta >  1.0) cosTheta =  1.0;
      if (cosTheta < -1.0) cosTheta = -1.0;
   }

   /*
   * Return value sin(theta) for precomputed cos(theta).
   */
   inline
   double Bend::sinTheta() const
   {  return sqrt(1.0 - cosTheta*cosTheta); }

   /*
   * Return theta in radians for precomputed cos(theta).
   */
   inline
   double Bend::theta() const
   {  return std::acos(cosTheta); }

} 
#endif
