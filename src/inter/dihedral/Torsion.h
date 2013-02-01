#ifndef INTER_TORSION_H
#define INTER_TORSION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>

#include <cmath>

namespace Inter
{

   using namespace Util;

   /**
   * Computes dihedral / torsion angle involving 3 bonds.
   *
   * Models the dihedral angle of 3 bond vectors b1, b2, and b3 
   * that connect 4 atoms. These three vectors must be passed to
   * the function computeAngle(), which computes the cosine of 
   * the resulting dihedral angle, and stores it in the public 
   * member cosPhi.
   *
   * See \ref Inter_Dihedral_Module in the file dihedral.mod for 
   * the definition of cosPhi.
   *
   * \ingroup Inter_Dihedral_Module
   */
   struct Torsion
   {
   
   public:

      /**
      * Cosine of dihedral angle.
      */
      double cosPhi;

      /**
      * Compute cosPhi.
      *
      * \param b1  bond vector from atom 0 to 1.
      * \param b2  bond vector from atom 1 to 2.
      * \param b3  bond vector from atom 2 to 3.
      */
      void computeAngle(const Vector& b1, const Vector& b2, const Vector& b3);
 
      /**
      * Return value of sin(phi) for precomputed cos(phi).
      */
      double sinPhi() const;

      /**
      * Return value of phi in radians for precomputed cos(phi).
      *
      * Returns principal value, in range 0 < phi < Pi.
      */
      double phi() const;

   };

   // Inline method definitions

   /* 
   * Calculate cosPhi.
   */ 
   inline void 
   Torsion::computeAngle(const Vector& b1, const Vector& b2, const Vector& b3)
   {
      Vector v1, v2;
      v1.cross(b1, b2);
      v2.cross(b2, b3);
      double d = sqrt(v1.square()*v2.square());
      cosPhi = v1.dot(v2)/d;
   }

   /*
   * Return value sin(phi) for precomputed cos(phi).
   */
   inline
   double Torsion::sinPhi() const
   {  return sqrt(1.0 - cosPhi*cosPhi); }

   /*
   * Return phi in radians for precomputed cos(phi).
   */
   inline
   double Torsion::phi() const
   {  return std::acos(cosPhi); }

} 
#endif
