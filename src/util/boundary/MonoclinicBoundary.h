#ifndef UTIL_MONOCLINIC_BOUNDARY_H
#define UTIL_MONOCLINIC_BOUNDARY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/crystal/LatticeSystem.h> // member
#include <util/containers/FSArray.h>    // member template
#include <util/space/Vector.h>          // member template argument
#include <util/space/IntVector.h>       // inline methods
#include <util/space/Dimension.h>       // member template argument
#include <util/global.h>                // asserts in inline methods

#include <cmath>
#include <iostream>

class MonoclinicBoundaryTest;

#ifndef UTIL_ORTHOGONAL
#define UTIL_ORTHOGONAL 0
#endif

namespace Util
{

   class Random;

   #ifdef UTIL_MPI
   template <class T> void send(MPI::Comm& comm, T& data, int dest, int tag);
   template <class T> void recv(MPI::Comm& comm, T& data, int source, int tag);
   template <class T> void bcast(MPI::Intracomm& comm, T& data, int root);
   #endif

   /**
   * A monoclinic periodic unit cell.
   *
   * \ingroup Boundary_Module
   */
   class MonoclinicBoundary 
   {

   public:

      /**
      * Constructor.
      */
      MonoclinicBoundary();

      /**
      * Set unit cell dimensions for orthorhombic boundary.
      *
      * Also sets all related lengths and volume.
      *
      * \param lengths  Vector of unit cell dimensions, orthogonal box.
      * \param t        displacement along z axis.
      */
      void setMonoclinic(const Vector &lengths, const double d);

      /**
      * Invalid function for monoclinic - throws Exception.
      *
      * \param lengths  Vector of unit cell dimensions, orthogonal box.
      */
      void setOrthorhombic(const Vector &lengths);

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
      * On output the vector r is shifted its image in the primary 
      * unit cell.
      *
      * \param r Vector of coordinates
      */
      void shift(Vector &r) const;

      /**
      * Shift Vector r to its image within the primary unit cell.
      *
      * This method shifts a vector to its image in the primary unit cell, 
      * and also increments a corresponding shift IntVector:
      *
      * If the Vector r is shifted by r -> r - \sum_i t[i]*a[i], then
      * the IntVector shift is shifted by  shift -> shift + t, where
      * t is an IntVector shift.
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
      * Generalized coordinates range from 0.0 <= Rg[i] < 1.0 within the
      * primitive cell, for i = 0,..,2.
      *
      * \param Rc Vector of Cartesian coordinates (input)
      * \param Rg Vector of generalized coordinates (output)
      */
      void transformCartToGen(const Vector& Rc, Vector& Rg) const;

      /**
      * Transform Vector of generalized coordinates to Cartesian coordinates.
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
      * Value is Monoclinic.
      */
      LatticeSystem latticeSystem();

      /**
      * Get Vector of distances between faces of primitive cell. 
      *
      * Each component of this vector is projection of a Bravais lattice unit
      * vector onto a line parallel to the corresponding reciprocal lattice
      * basis vector.  The resulting distances are distances between faces
      * faces of the primitive unit cell, which are normal to the reciprocal
      * lattice unit vectors. 
      */
      const Vector& lengths() const;

      /**
      * Get distance across primitive cell parallel to reciprocal axis i.
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

      /// Minimum coordinates for orthogonal unit cell
      Vector minima_;

      /// Maximum coordinates for orthogonal unit cell
      Vector maxima_;

      /// Lengthsof orthogonal unit cell:  l_[i] = maxima_[i] - minima_[i].
      Vector l_;

      /// Displacement (tilt) of primitive cell in the z-direction.
      double tilt_;

      /// Half region lengths: halfL_[i] = 0.5*l_[i].
      Vector halfL_;

      /// Half region lengths: invL_[i] = 1.0/l_[i].
      Vector invL_;

      /// Lengths of the primitive projected along reciprocal basis vectors.
      Vector lengths_;

      /// Volume: V = l_[0]*l_[1]*l_[2].
      double volume_;

      /// constants used for distancing: u2 = dz + c3_*dy.
      double c3_;

      /// Length of bravais basis vector 1 (tilted).
      double e_;

      /// Minimum distance across the unit cell.
      double minLength_;

      /**
      * Array of Bravais lattice vectors.
      */
      FSArray<Vector, Dimension>  bravaisBasisVectors_;

      /**
      * Array of Reciprocal lattice vectors.
      */
      FSArray<Vector, Dimension>  reciprocalBasisVectors_;

      /**
      * Actual lattice system (Monoclinic)
      */
      LatticeSystem lattice_;

   // friends:

      /// Unit test
      friend class ::MonoclinicBoundaryTest;

      /// istream extractor
      friend std::istream& operator >> (std::istream& in, 
                                        MonoclinicBoundary& boundary);

      /// ostream inserter
      friend std::ostream& operator << (std::ostream& out, 
                                        const MonoclinicBoundary& boundary);

      #ifdef UTIL_MPI
      friend void send<>(MPI::Comm& comm, MonoclinicBoundary& data, 
                         int dest, int tag);

      friend void recv<>(MPI::Comm& comm, MonoclinicBoundary& data, 
                         int source, int tag);

      friend void bcast<>(MPI::Intracomm& comm, MonoclinicBoundary& data, 
                          int root);
      #endif

      /**
      * Reset all quantities that depend upon lengths.
      */ 
      void reset();

   };

   // Inline method definitions

   /* 
   * Return Vector of lengths by const reference.
   */
   inline const Vector& MonoclinicBoundary::lengths() const 
   {  return lengths_; }

   /* 
   * Get length = maximum - minimum in direction i.
   */
   inline double MonoclinicBoundary::length(int i) const 
   {  return lengths_[i]; }

   /* 
   * Return region volume.
   */
   inline double MonoclinicBoundary::volume() const
   {  return volume_; }

   /* 
   * Return Bravais lattice basis vector number i.
   */
   inline const Vector& MonoclinicBoundary::bravaisBasisVector(int i) const
   {  return bravaisBasisVectors_[i]; }

   /* 
   * Return reciprocal lattice basis vector number i.
   */
   inline const Vector& MonoclinicBoundary::reciprocalBasisVector(int i) const
   {  return reciprocalBasisVectors_[i]; }

   /* 
   * Return the maximum validity range of the distances.
   */
   inline double MonoclinicBoundary::minLength() const
   {  return minLength_; }

   /* 
   * Shift Vector r to periodic cell, R[axis] < r[axis] < maxima_[axis].
   */
   inline void MonoclinicBoundary::shift(Vector& r) const
   {
       if(r[0] >= maxima_[0]) {
          r[0] -= l_[0];
          assert(r[0] < maxima_[0]);
       } else
       if (r[0] <  minima_[0]) {
          r[0] += l_[0];
          assert(r[0] >= minima_[0]);
       }
       if(r[1] >= maxima_[1]) {
          r[1] -= l_[1];
          r[2] -= tilt_;
          assert(r[1] < l_[1]);
       } else
       if (r[1] <  minima_[1]) {
          r[1] += l_[1];
          r[2] += tilt_;
          assert(r[1] >= minima_[1]);
       }
       double u2 = r[2] + c3_*r[1];
       if(u2 >= maxima_[2]) {
          r[2] -= l_[2];
          assert(r[2] + c3_*r[1] >= minima_[2]);
       } else
       if (u2 <  minima_[2]) {
          r[2] += l_[2];
          assert(r[2] + c3_*r[1] >= minima_[2]);
       }
   }

   /* 
   * Shift Vector r to periodic cell, R[axis] < r[axis] < maxima_[axis].
   */
   inline void MonoclinicBoundary::shift(Vector& r, IntVector& shift) const
   {
       if(r[0] >= maxima_[0]) {
          r[0] -= l_[0];
          ++(shift[0]);
          assert(r[0] < maxima_[0]);
       } else
       if (r[0] <  minima_[0]) {
          r[0] += l_[0];
          --(shift[0]);
          assert(r[0] >= minima_[0]);
       }
       if (r[1] >= maxima_[1]) {
          r[1] -= l_[1];
          r[2] -= tilt_;
          ++(shift[1]);
          assert(r[1] < l_[1]);
       } else
       if (r[1] <  minima_[1]) {
          r[1] += l_[1];
          r[2] += tilt_;
          --(shift[1]);
          assert(r[1] >= minima_[1]);
       }
       double u2 = r[2] + c3_*r[1];
       if (u2 >= maxima_[2]) {
          r[2] -= l_[2];
          ++(shift[2]);
          assert(r[2] + c3_*r[1] < maxima_[2]);
       } else
       if (u2 <  minima_[2]) {
          r[2] += l_[2];
          --(shift[2]);
          assert(r[2] + c3_*r[1] >= minima_[2]);
       }
   }

   /* 
   * Calculate squared distance by half-parallelogram convention.
   */
   inline 
   double MonoclinicBoundary::distanceSq(const Vector &r1, const Vector &r2, 
                               IntVector& translate) const
   {

      double dx = r1[0] - r2[0];
      if ( fabs(dx) > halfL_[0] ) {
         if (dx > 0.0) {
            dx -= l_[0];
            translate[0] = -1;
         } else {
            dx += l_[0];
            translate[0] = +1;
         }
         assert(fabs(dx) <= halfL_[0]);
      } else {
         translate[0] = 0;
      }

      double dy = r1[1] - r2[1];
      double dz = r1[2] - r2[2];
      double u2 = dz + c3_*dy;
      if ( fabs(dy) > halfL_[1] ) {
         if (dy > 0.0) {
            dy -= l_[1];
            dz -= tilt_;
            translate[1] = -1;
         } else {
            dy += l_[1];
            dz += tilt_;
            translate[1] = +1;
         }
         assert(fabs(dy) <= halfL_[1]);
      } else {
         translate[1] = 0;
      }

      if ( fabs(u2) >  halfL_[2]) {
         if (u2 > 0.0) {
            dz -= l_[2];
            translate[2] = -1;
         } else {
            dz += l_[2];
            translate[2] = +1;
         }
         assert(fabs(dz + c3_*dy) <= halfL_[2]);
      } else {
         translate[2] = 0;
      }
      return dx*dx+dy*dy+dz*dz;
   }

   /* 
   * Return squared distance and separation by half-parallelogram convention.
   */
   inline 
   double MonoclinicBoundary::distanceSq(const Vector &r1, 
                                         const Vector &r2) const
   {

      double dx = r1[0] - r2[0];
      if ( fabs(dx) > halfL_[0] ) {
         if (dx > 0.0) {
            dx -= l_[0];
         } else {
            dx += l_[0];
         }
         assert(fabs(dx) <= halfL_[0]);
      }

      double dy = r1[1] - r2[1];
      double dz = r1[2] - r2[2];
      double u2 = dz + c3_*dy;
      if ( fabs(dy) > halfL_[1] ) {
         if (dy > 0.0) {
            dy -= l_[1];
            dz -= tilt_;
         } else {
            dy += l_[1];
            dz += tilt_;
         }
         assert(fabs(dy) <= halfL_[1]);
      }

      if ( fabs(u2) >  halfL_[2]) {
         if (u2 > 0.0) {
            dz -= l_[2];
         } else {
            dz += l_[2];
         }
         assert(fabs(dz+c3_*dy) <= halfL_[2]);
      }

      return dx*dx+dy*dy+dz*dz;
   }

   /* 
   * Calculate squared distance between two positions, and separation vector,
   * by half-parallelogram convention.
   */
   inline 
   double MonoclinicBoundary::distanceSq(const Vector &r1, const Vector &r2, 
                                         Vector &dr) const
   {

      double dx = r1[0] - r2[0];
      if (fabs(dx) > halfL_[0]) {
         if (dx > 0.0) {
            dx -= l_[0];
         } else {
            dx += l_[0];
         }
         assert(fabs(dx) <= halfL_[0]);
      }

      double dy = r1[1] - r2[1];
      double dz = r1[2] - r2[2];
      double u2 = dz + c3_*dy;

      if (fabs(dy) > halfL_[1]) {
         if (dy > 0.0) {
            dy -= l_[1];
            dz -= tilt_;
         } else {
            dy += l_[1];
            dz += tilt_;
         }
         assert(fabs(dy) <= halfL_[1]);
      }

      if (fabs(u2) >  halfL_[2]) {
         if (u2 > 0.0) {
            dz -= l_[2];
         } else {
	    dz += l_[2];
         }
         assert(fabs(dz+c3_*dy) <= halfL_[2]);
      }

      dr[0] = dx;
      dr[1] = dy;
      dr[2] = dz;

      return (dx*dx + dy*dy + dz*dz);
   }

   /**
   * istream extractor for a MonoclinicBoundary.
   *
   * \param  in        input stream
   * \param  boundary  MonoclinicBoundary to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, MonoclinicBoundary& boundary);

   /**
   * ostream inserter for an MonoclinicBoundary.
   *
   * \param  out      output stream
   * \param  boundary MonoclinicBoundary to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, 
                              const MonoclinicBoundary& boundary);

   /**
   * Return actual lattice system.
   */
   inline LatticeSystem MonoclinicBoundary::latticeSystem()
   { return lattice_; }

   inline 
   void MonoclinicBoundary::transformCartToGen(const Vector& Rc, Vector& Rg) 
   const
   {
      Rg[0] = Rc[0] * invL_[0];
      Rg[1] = Rc[1] * invL_[1];
      Rg[2] = (Rc[2] + c3_ *Rc[1]) * invL_[2];
   }
      
   /**
   * Returns the Generalized coordinates of Rc in Rg.
   */
   inline 
   void MonoclinicBoundary::transformGenToCart(const Vector& Rg, Vector& Rc) 
   const
   {
      Rc[0] = Rg[0]*l_[0];
      Rc[1] = Rg[1]*l_[1];
      Rc[2] = Rg[2]*l_[2] + Rg[1]*tilt_;
   }

}
 
#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiTraits.h>

namespace Util
{

   /**
   * Send an MonoclinicBoundary via MPI.
   */
   template <>
   void send<MonoclinicBoundary>(MPI::Comm& comm, MonoclinicBoundary& data, 
                                 int dest, int tag);

   /**
   * Receive an MonoclinicBoundary via MPI.
   */
   template <>
   void recv<MonoclinicBoundary>(MPI::Comm& comm, MonoclinicBoundary& data, 
                                 int source, int tag);

   /**
   * Broadcast an MonoclinicBoundary via MPI.
   */
   template <>
   void bcast<MonoclinicBoundary>(MPI::Intracomm& comm, 
                                  MonoclinicBoundary& data, int root);

   /**
   * Explicit specialization MpiTraits<MonoclinicBoundary>.
   */
   template <>
   class MpiTraits<MonoclinicBoundary>
   {
   public:
      static MPI::Datatype type;         ///< MPI Datatype
      static bool hasType;               ///< Is the MPI type initialized?
   };

}
#endif // ifdef  UTIL_MPI
#endif // ifndef UTIL_MONOCLINIC_BOUNDARY_H
