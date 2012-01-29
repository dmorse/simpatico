#ifndef ORTHO_BOUNDARY_BASE_H
#define ORTHO_BOUNDARY_BASE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthoRegion.h"              // base class
#include <util/containers/FSArray.h>  // member
#include <util/space/Vector.h>        // member 
#include <util/space/Dimension.h>     // member template parameter
#include <util/space/IntVector.h>     // inline methods
#include <util/global.h>              // asserts in inline methods

#include <cmath>

class OrthoBoundaryTest;
namespace Util{ class Random; }

namespace DdMd
{

   using namespace Util;

   /**
   * Base class for othorhombic periodic simulation cells.
   *
   * This class is designed as a base class for classes that
   * represent Orthorhombic, Tetragonal, and Cubic periodic
   * unit cells. It provides methods to impose periodic boundary 
   * conditions, and to calculate separations and distances
   * using the minimum image convention, for a general
   * orthorhombic unit cell. 
   *
   * \ingroup Boundary_Module
   */
   class OrthoBoundaryBase : protected OrthoRegion
   {
   
   public:

      /**
      * Constructor.
      */
      OrthoBoundaryBase();

      ///\name Accessors
      //@{

      /**
      * Get Vector of lengths by const reference.
      */
      const Vector& lengths() const;

      /**
      * Get length in Cartesian direction i.
      *
      * \param i index of Cartesian direction, 0 <= i < Dimension
      */
      double length(int i) const;

      /**
      * Return unit cell volume.
      */
      double volume() const;

      /**
      * Return Bravais lattice vector i
      *  
      * \param i basis Vector index.
      */ 
      const Vector& bravaisBasisVector(int i) const;

      /**
      * Return reciprocal lattice basis vector i
      *  
      * \param i basis Vector index.
      */ 
      const Vector& reciprocalBasisVector(int i) const;

      /**
      * Generate random position within the primary unit cell.
      *
      * Generates coordinates r[i], i=1,..3 with minima_[i] < r[i] < maxima_[i].
      *
      * \param random random number generator object
      * \param r      Vector of 3 random coordinates (upon return)
      */
      void randomPosition(Random &random, Vector &r) const;

      //@}
      ///\name Periodic Boundary Conditions
      //@{

      /**
      * Shift Vector r to its image within the primary unit cell.
      *
      * One output, each coordinate r[i] is shifted by a multiple of length[i]
      * so as to lie within the range minima_[i] < r[i] < maxima_[i].
      *
      * Precondition: The algorithm assumes that on input, for each i=0,..,2,
      * minima_[i] - lengths_[i] < r[i] < maxima_[i] + lengths_[i]
      *
      * \param r Vector of coordinates
      */
      void shift(Vector &r) const;

      #if 0
      /**
      * Shift Vector r to its image within the primary unit cell.
      *
      * This metho has the same effect on r as shift(Vector& r). On
      * return, the elements of the extra parameter translate indicate
      * the translation required to bring r into the primary unit cell:
      * On return, r has components r[i] -> r[i] + translate[i]*length[i].
      *
      * \param r         Vector of coordinates
      * \param translate integers multiplied by boundary lengths
      */
      void shift(Vector &r, IntVector& translate) const;
      #endif

      /**
      * Return square distance between positions r1 and r2, using the nearest
      * the nearest image convention for the separation Vector.
      *
      * \param r1  first position Vector
      * \param r2  second position Vector
      * \return square of distance between r1 and r2, using nearest image.
      */
      double distanceSq(const Vector &r1, const Vector &r2) const;

      /**
      * Return square distance between positions r1 and r2, using the nearest
      * the nearest image convention for the separation Vector.
      *
      * \param r1    first position Vector
      * \param r2    second position Vector
      * \param shift shift added to r1 to create nearest image of r2.
      * \return square of distance between r1 and r2, using nearest image.
      */
      double distanceSq(const Vector &r1, const Vector &r2, IntVector& shift) 
      const;

      /**
      * Return the squared distance between positions r1 and r2, using the
      * nearest image convention, and calculate the separation Vector.
      *
      * Upon return, Vector dr contains the separation r1 - r2, using
      * the nearest image convention. Returns the square of the absolute
      * magnitude of the separation dr.
      *
      * \param r1  first position Vector
      * \param r2  second position Vector
      * \param dr  separation Vector (upon return)
      * \return square of separation Vector dr
      */
      double distanceSq(const Vector &r1, const Vector &r2, Vector &dr) const;

      //@}

   protected:

      /**
      * Reset all quantities that depend upon lengths.
      */ 
      void reset();

   private:

      FSArray<Vector, Dimension>  bravaisBasisVectors_;

      FSArray<Vector, Dimension>  reciprocalBasisVectors_;

   // friends:

      /// Unit test
      friend class ::OrthoBoundaryTest;

   };

   // Inline method definitions

   /* 
   * Return Vector of lengths by const reference.
   */
   inline const Vector& OrthoBoundaryBase::lengths() const 
   {  return lengths_; }

   /* 
   * Get length = maximum - minimum in direction i.
   */
   inline double OrthoBoundaryBase::length(int i) const 
   {  return lengths_[i]; }

   /* 
   * Return region volume.
   */
   inline double OrthoBoundaryBase::volume() const
   {  return volume_; }

   /* 
   * Return Bravais lattice basis vector number i.
   */
   inline const Vector& OrthoBoundaryBase::bravaisBasisVector(int i) const
   {  return bravaisBasisVectors_[i]; }

   /* 
   * Return reciprocal lattice basis vector number i.
   */
   inline const Vector& OrthoBoundaryBase::reciprocalBasisVector(int i) const
   {  return reciprocalBasisVectors_[i]; }

   /* 
   * Shift Vector r to periodic cell, R[axis] < r[axis] < maxima_[axis].
   */
   inline void OrthoBoundaryBase::shift(Vector& r) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if( r[i] >= maxima_[i] ) {
            r[i] = r[i] - lengths_[i];
         } else
         if ( r[i] <  minima_[i] ) {
            r[i] = r[i] + lengths_[i];
         }
      }
   }

   #if 0
   /* 
   * Shift Vector r to periodic cell, R[axis] < r[axis] < maxima_[axis].
   */
   inline void OrthoBoundaryBase::shift(Vector& r, IntVector& translate) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if ( r[i] >= maxima_[i] ) {
            r[i] = r[i] - lengths_[i];
            translate[i] = -1;
         } else
         if ( r[i] <  minima_[i] ) {
            r[i] = r[i] + lengths_[i];
            translate[i] = +1;
         } else {
            translate[i] = 0;
         }
      }
   }
   #endif

   /* 
   * Calculate squared distance by minimum image convention.
   */
   inline 
   double OrthoBoundaryBase::distanceSq(const Vector &r1, const Vector &r2, 
                               IntVector& translate) const
   {
      double dr;
      double norm = 0.0;
      for (int i = 0; i < Dimension; ++i) {
         dr = r1[i] - r2[i];
         if ( fabs(dr) > halfLengths_[i] ) {
            if (dr > 0.0) {
               dr -= lengths_[i];
               translate[i] = -1;
            } else {
               dr += lengths_[i];
               translate[i] = +1;
            }
            assert(fabs(dr) <= halfLengths_[i]);
         } else {
            translate[i] = 0;
         }
         norm += dr*dr;
      }
      return norm;
   }

   /* 
   * Return squared distance and separation with minimum image convention.
   */
   inline 
   double OrthoBoundaryBase::distanceSq(const Vector &r1, const Vector &r2) const
   {
      double dr;
      double norm = 0.0;
      for (int i = 0; i < Dimension; ++i) {
         dr = r1[i] - r2[i];
         if ( fabs(dr) > halfLengths_[i] ) {
            if ( dr > 0.0 ) {
               dr -= lengths_[i];
            } else {
               dr += lengths_[i];
            }
            assert(fabs(dr) <= halfLengths_[i]);
         }
         norm += dr*dr;
      }
      return norm;
   }

   /* 
   * Calculate squared distance between two positions, and separation vector,
   * using the minimum image convention.
   */
   inline 
   double OrthoBoundaryBase::distanceSq(const Vector &r1, const Vector &r2, 
                                          Vector &dr) const
   {
      for (int i = 0; i < Dimension; ++i) {
         dr[i] = r1[i] - r2[i];
         if (fabs(dr[i]) > halfLengths_[i]) {
            if (dr[i] > 0.0) {
               dr[i] -= lengths_[i];
            } else {
               dr[i] += lengths_[i];
            }
            assert(fabs(dr[i]) <= halfLengths_[i]);
         }
      }
      return ( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );
   }

} 
#endif
