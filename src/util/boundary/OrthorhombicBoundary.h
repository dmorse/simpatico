#ifndef ORTHORHOMBIC_BOUNDARY_H
#define ORTHORHOMBIC_BOUNDARY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthoRegion.h"                // base class
#include <util/crystal/LatticeSystem.h> // member
#include <util/containers/FSArray.h>    // member template
#include <util/space/Vector.h>          // member template argument
#include <util/space/IntVector.h>       // inline methods
#include <util/space/Dimension.h>       // member template argument
#include <util/global.h>                // asserts in inline methods

#include <cmath>
#include <iostream>

class OrthorhombicBoundaryTest;

#ifndef UTIL_ORTHOGONAL
#define UTIL_ORTHOGONAL 1
#endif

namespace Util
{

   class Random;

   /**
   * An orthorhombic periodic unit cell.
   *
   * \ingroup Boundary_Module
   */
   class OrthorhombicBoundary : private OrthoRegion
   {

   public:

      /**
      * Constructor.
      */
      OrthorhombicBoundary();

      /**
      * Set unit cell dimensions for orthorhombic boundary.
      *
      * Also sets all related lengths and volume.
      *
      * \param lengths  Vector of unit cell lengths
      */
      void setLengths(const Vector &lengths);

      /**
      * Set unit cell dimensions for tetragonal boundary.
      *
      * \param ab unit cell dimensions in x and y directions
      * \param c  unit cell length in z direction
      */
      void setTetragonalLengths(double ab, double c);

      /**
      * Set unit cell dimensions for a cubic boundary.
      *
      * \param a unit cell length in x, y, and z directions.
      */
      void setCubicLengths(double a);

      /**
      * Serialize to/from an archive.
      * 
      * \param ar       saving or loading archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

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

      /**
      * Shift Vector r to its image within the primary unit cell.
      *
      * This method maps an atomic position to lie in the primary cell, 
      * and also increments the atomic shift IntVector:
      *
      * If r[i] -> r[i] - t*length_[i], then shift[i] -> shift[i] + t.
      *
      * \sa Atom:shift()
      *
      * \param r     Vector of coordinates
      * \param shift integer shifts required to obtain "true" coordinates.
      */
      void shift(Vector &r, IntVector& shift) const;

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
      ///\name Coordinate Transformations
      //@{

      /**
      * Transform Cartesian Vector to generalized coordinates.
      *
      * Generalized coordinates range from 0.0 < Rg[i] < 1.0 within the
      * primitive cell, for i=0,..,2.
      *
      * \param Rc Vector of Cartesian coordinates (input)
      * \param Rg Vector of generalized coordinates (output)
      */
      void transformCartToGen(const Vector& Rc, Vector& Rg) const;

      /**
      * Transform Vector of generalized coordinates to Cartesian Vector.
      *
      * \param Rg Vector of generalized coordinates (input)
      * \param Rc Vector of Cartesian coordinates (output)
      */
      void transformGenToCart(const Vector& Rg, Vector& Rc) const;

      //@}
      ///\name Accessors
      //@{

      /**
      * Return actual lattice system.
      *
      * Value can be Cubic, Tetragonal, or Orthorhombic.
      */
      LatticeSystem latticeSystem();

      /**
      * Get Vector of unit cell lengths by const reference.
      */
      const Vector& lengths() const;

      /**
      * Get length in Cartesian direction i.
      *
      * \param i index of Cartesian direction, 0 <= i < Dimension
      */
      double length(int i) const;

      /**
      * Get minimum length across primitive unit cell.
      */
      double minLength() const;

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
      * Generates Vector r[i], i=1,..,3 with minima_[i] < r[i] < maxima_[i].
      *
      * \param random random number generator object
      * \param r      Vector of random coordinates (upon return)
      */
      void randomPosition(Random &random, Vector &r) const;

      /**
      * Return true if valid, or throw Exception.
      */
      bool isValid();

      //@}

   private:

      /**
      * Array of Bravais lattice vectors.
      */
      FSArray<Vector, Dimension>  bravaisBasisVectors_;

      /**
      * Array of Reciprocal lattice vectors.
      */
      FSArray<Vector, Dimension>  reciprocalBasisVectors_;

      /**
      * Vector of inverse box dimensions.
      */
      Vector invLengths_;

      /**
      * Minimum distance across the unit cell.
      */
      double minLength_;

      /**
      * Actual lattice system (Orthorhombic, Tetragonal, or Cubic)
      */
      LatticeSystem lattice_;

   // friends:

      /// Unit test
      friend class ::OrthorhombicBoundaryTest;

      /// istream extractor
      friend std::istream& operator >> (std::istream& in, 
                                        OrthorhombicBoundary& boundary);

      /// ostream inserter
      friend std::ostream& operator << (std::ostream& out, 
                                        const OrthorhombicBoundary& boundary);

      /**
      * Reset all quantities that depend upon lengths.
      */ 
      void reset();

   };

   // Inline method definitions

   /* 
   * Return Vector of lengths by const reference.
   */
   inline const Vector& OrthorhombicBoundary::lengths() const 
   {  return lengths_; }

   /* 
   * Get length = maximum - minimum in direction i.
   */
   inline double OrthorhombicBoundary::length(int i) const 
   {  return lengths_[i]; }

   /* 
   * Return the maximum validity range of the distances.
   */
   inline double OrthorhombicBoundary::minLength() const
   {  return minLength_; }

   /* 
   * Return region volume.
   */
   inline double OrthorhombicBoundary::volume() const
   {  return volume_; }

   /* 
   * Return Bravais lattice basis vector number i.
   */
   inline const Vector& OrthorhombicBoundary::bravaisBasisVector(int i) const
   {  return bravaisBasisVectors_[i]; }

   /* 
   * Return reciprocal lattice basis vector number i.
   */
   inline const Vector& OrthorhombicBoundary::reciprocalBasisVector(int i) const
   {  return reciprocalBasisVectors_[i]; }

   /* 
   * Shift Vector r to periodic cell, R[axis] < r[axis] < maxima_[axis].
   */
   inline void OrthorhombicBoundary::shift(Vector& r) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if( r[i] >= maxima_[i] ) {
            r[i] = r[i] - lengths_[i];
            assert(r[i] < maxima_[i]);
         } else
         if ( r[i] <  minima_[i] ) {
            r[i] = r[i] + lengths_[i];
            assert(r[i] >= minima_[i]);
         }
      }
   }

   /* 
   * Shift Vector r to periodic cell, R[axis] < r[axis] < maxima_[axis].
   */
   inline void OrthorhombicBoundary::shift(Vector& r, IntVector& shift) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if (r[i] >= maxima_[i]) {
            r[i] = r[i] - lengths_[i];
            ++(shift[i]);
            assert(r[i] < maxima_[i]);
         } else
         if (r[i] <  minima_[i]) {
            r[i] = r[i] + lengths_[i];
            --(shift[i]);
            assert(r[i] >= minima_[i]);
         }
      }
   }

   /* 
   * Calculate squared distance by minimum image convention.
   */
   inline 
   double OrthorhombicBoundary::distanceSq(const Vector &r1, const Vector &r2, 
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
   double OrthorhombicBoundary::distanceSq(const Vector &r1, const Vector &r2) const
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
   double OrthorhombicBoundary::distanceSq(const Vector &r1, const Vector &r2, 
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

   /**
   * istream extractor for a OrthorhombicBoundary.
   *
   * \param  in        input stream
   * \param  boundary  OrthorhombicBoundary to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, OrthorhombicBoundary& boundary);

   /**
   * ostream inserter for an OrthorhombicBoundary.
   *
   * \param  out      output stream
   * \param  boundary OrthorhombicBoundary to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, 
                              const OrthorhombicBoundary& boundary);

   /**
   * Return actual lattice system.
   */
   inline LatticeSystem OrthorhombicBoundary::latticeSystem()
   { return lattice_; }

   /*
   * Transform Cartesian Vector Rc to generalized Vector Rg.
   */
   inline void 
   OrthorhombicBoundary::transformCartToGen(const Vector& Rc, Vector& Rg) 
   const
   {
      for (int i = 0; i < Dimension; ++i) {
         Rg[i] = Rc[i] * invLengths_[i];
      }
   }
      
   /*
   * Transform Cartesian Vector Rc to generalized Vector Rg.
   */
   inline void 
   OrthorhombicBoundary::transformGenToCart(const Vector& Rg, Vector& Rc) 
   const
   {
      for (int i = 0; i < Dimension; ++i) {
         Rc[i] = Rg[i] * lengths_[i];
      }
   }

}
 
#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiTraits.h>

namespace Util
{

   /**
   * Send an OrthorhombicBoundary via MPI.
   */
   template <>
   void send<Util::OrthorhombicBoundary>(MPI::Comm& comm, 
             Util::OrthorhombicBoundary& data, int dest, int tag);

   /**
   * Receive an OrthorhombicBoundary via MPI.
   */
   template <>
   void recv<Util::OrthorhombicBoundary>(MPI::Comm& comm, 
             Util::OrthorhombicBoundary& data, int source, int tag);

   /**
   * Broadcast an OrthorhombicBoundary via MPI.
   */
   template <>
   void bcast<Util::OrthorhombicBoundary>(MPI::Intracomm& comm, 
              Util::OrthorhombicBoundary& data, int root);

   /**
   * Explicit specialization MpiTraits<OrthorhombicBoundary>.
   */
   template <>
   class MpiTraits<Util::OrthorhombicBoundary>
   {
   public:
      static MPI::Datatype type;         ///< MPI Datatype
      static bool hasType;               ///< Is the MPI type initialized?
   };

}
#endif // ifdef  UTIL_MPI

#endif // ifndef ORTHORHOMBIC_BOUNDARY_H
