#ifndef MONOCLINIC_BOUNDARYMI_H
#define MONOCLINIC_BOUNDARYMI_H

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

class MonoclinicBoundaryMITest;

namespace Util
{

   class Random;

   /**
   * An orthorhombic periodic unit cell.
   *
   * \ingroup Boundary_Module
   */
   class MonoclinicBoundaryMI 
   {

   public:

      /**
      * Constructor.
      */
      MonoclinicBoundaryMI();

      /**
      * Set unit cell dimensions for orthorhombic boundary.
      *
      * Also sets all related lengths and volume.
      *
      * \param lengths  Vector of unit cell lengths
      */
      void setLengths(const Vector &lengths, const double d);

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
      ///\name Accessors
      //@{

      /**
      * Return actual lattice system.
      *
      * Value is Monoclinic.
      */
      LatticeSystem latticeSystem();

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

      /// Minimum coordinates: Require r[i] >= minima_[i].
      Vector minima_;

      /// Maximum coordinates: Require r[i] <  maxima_[i].
      Vector maxima_;

      /// Box lengths:  lengths_[i] = maxima_[i] - minima_[i].
      Vector lengths_;

      /// Box tilt_ in Monoclinic Box. tilt_ is in the z-direction.
      double tilt_;

      /// Half region lengths: halfLengths_[i] = 0.5*lengths_[i].
      Vector halfLengths_;

      /// Volume: V = lengths_[0]*lengths_[1]*lengths_[2].
      double volume_;

      /// constants used for distancing: u1=c1*dy.
      double c1_;

      /// constants used for distancing: u2=c2*dz+c3*dy.
      double c2_;

      /// constants used for distancing: u2=c2*dz+c3*dy.
      double c3_;

      /// the length of the minor axis of the Monoclinic parallelogram.
      double e_;

      /// the half length of the minor axis of the Monoclinic parallelogram.
      double halfe_;

      /// the length of the minor diameter of the Monoclinic parallelogram.
      double d1_;

      /// the half length of the minor diameter of the Monoclinic parallelogram.
      double halfd1_;

      /// the length of the major diameter of the Monoclinic parallelogram.
      double d2_;

      /// the half length of major diameter of the Monoclinic parallelogram.
      double halfd2_;

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
      friend class ::MonoclinicBoundaryMITest;

      /// istream extractor
      friend std::istream& operator >> (std::istream& in, 
                                        MonoclinicBoundaryMI& boundary);

      /// ostream inserter
      friend std::ostream& operator << (std::ostream& out, 
                                        const MonoclinicBoundaryMI& boundary);

      /**
      * Reset all quantities that depend upon lengths.
      */ 
      void reset();

   };

   // Inline method definitions

   /* 
   * Return Vector of lengths by const reference.
   */
   inline const Vector& MonoclinicBoundaryMI::lengths() const 
   {  return lengths_; }

   /* 
   * Get length = maximum - minimum in direction i.
   */
   inline double MonoclinicBoundaryMI::length(int i) const 
   {  return lengths_[i]; }

   /* 
   * Return region volume.
   */
   inline double MonoclinicBoundaryMI::volume() const
   {  return volume_; }

   /* 
   * Return Bravais lattice basis vector number i.
   */
   inline const Vector& MonoclinicBoundaryMI::bravaisBasisVector(int i) const
   {  return bravaisBasisVectors_[i]; }

   /* 
   * Return reciprocal lattice basis vector number i.
   */
   inline const Vector& MonoclinicBoundaryMI::reciprocalBasisVector(int i) const
   {  return reciprocalBasisVectors_[i]; }

   /* 
   * Shift Vector r to periodic cell, R[axis] < r[axis] < maxima_[axis].
   */
   inline void MonoclinicBoundaryMI::shift(Vector& r) const
   {
         if( r[0] >= maxima_[0] ) {
            r[0] = r[0] - lengths_[0];
            assert(r[0] < maxima_[0]);
         } else
         if ( r[0] <  minima_[0] ) {
            r[0] = r[0] + lengths_[0];
            assert(r[0] >= minima_[0]);
         }
         if( r[1] >= maxima_[1] ) {
            r[1] = r[1] - lengths_[1];
	    r[2] = r[2] - tilt_;
            assert(r[1] < maxima_[1]);
         } else
         if ( r[1] <  minima_[1] ) {
            r[1] = r[1] + lengths_[1];
	    r[2] = r[2] - tilt_;
            assert(r[1] >= minima_[1]);
         }
         if( r[2] >= maxima_[2] ) {
            r[2] = r[2] - lengths_[2];
            assert(r[2] < maxima_[2]);
         } else
         if ( r[2] <  minima_[2] ) {
            r[2] = r[2] + lengths_[2];
            assert(r[2] >= minima_[2]);
         }
   }

   /* 
   * Shift Vector r to periodic cell, R[axis] < r[axis] < maxima_[axis].
   */
   inline void MonoclinicBoundaryMI::shift(Vector& r, IntVector& shift) const
   {
         if( r[0] >= maxima_[0] ) {
            r[0] = r[0] - lengths_[0];
	    ++(shift[0]);
            assert(r[0] < maxima_[0]);
         } else
         if ( r[0] <  minima_[0] ) {
            r[0] = r[0] + lengths_[0];
	    --(shift[0]);
            assert(r[0] >= minima_[0]);
         }
         if( r[1] >= maxima_[1] ) {
            r[1] = r[1] - lengths_[1];
	    r[2] = r[2] - tilt_;
	    ++(shift[1]);
            assert(r[1] < maxima_[1]);
         } else
         if ( r[1] <  minima_[1] ) {
            r[1] = r[1] + lengths_[1];
	    r[2] = r[2] - tilt_;
	    --(shift[1]);
            assert(r[1] >= minima_[1]);
         }
         if( r[2] >= maxima_[2] ) {
            r[2] = r[2] - lengths_[2];
            assert(r[2] < maxima_[2]);
	    ++(shift[2]);
         } else
         if ( r[2] <  minima_[2] ) {
            r[2] = r[2] + lengths_[2];
            assert(r[2] >= minima_[2]);
	    --(shift[2]);
         }
   }

   inline 
   double MonoclinicBoundaryMI::distanceSq(const Vector &r1, const Vector &r2, 
                               IntVector& translate) const
   {
      double dx = r1[0] - r2[0];
      double dy = r1[1] - r2[1];
      double dz = r1[2] - r2[2];
      double norm = 0.0;
      translate[0] = 0;
      translate[1] = 0;
      translate[2] = 0;

      if ( fabs(dx) > halfLengths_[0] ) {
         if (dx > 0.0) {
            dx -= lengths_[0];
            translate[0] += -1;
         } else {
            dx += lengths_[0];
            translate[0] += +1;
         }
         assert(fabs(dx) <= halfLengths_[0]);
      }

      bool done = true;
      double h1 = 0.0;
      double h2 = 0.0;
      double h3 = 0.0; 
      
      while(done){

	 h1 = (dy*lengths_[1]+dz*tilt_) / e_;
         if ( fabs(h1) > halfe_ ) {
	    done = true;
            if (h1 > 0.0) {
               dy -= lengths_[1];
	       dz -= tilt_;
               translate[1] += -1;
            } else {
               dy += lengths_[1];
	       dz += tilt_;
               translate[1] += +1;
            }
            assert(fabs(dy*length_[1]+dz*tilt) / e_ <= halfe_);
         } else {
	    done = false;
         }
         
         if ( fabs(dz) > halfLengths_[2] ) {
	    done = true;
            if (dz > 0.0) {
               dz -= lengths_[2];
               translate[2] += -1;
            } else {
               dz += lengths_[2];
               translate[0] += +1;
            }
            assert(fabs(dz) <= halfLengths_[2]);
         } else {
	    done = false;
         }

	 if ( (tilt_ > 0.0) ) {
         h2 = (dy*tilt_+dz*(tilt_-lengths_[2])) / d1_; 
         if ( fabs(h2) > halfd1_ ) {
	    done = true;
            if (h2 > 0.0) {
	       dy -= lengths_[1];
               dz -= tilt_-lengths_[2];
	       translate[1] += -1;
               translate[2] += +1;
            } else {
	       dy += lengths_[1];
               dz += tilt_-lengths_[2];
	       translate[1] += +1;
               translate[2] += -1;
            }
	    h2 = (dy*tilt_+dz*(tilt_-lengths_[2])) / d1_;
            assert(fabs(h2) <= halfd1_);
         } else {
	    done = false;
         }
         } else {
         h3 = (dy*tilt_+dz*(tilt_+lengths_[2])) / d2_;
         if ( fabs(h3) > halfd2_ ) {
	    done = true;
            if (h3 > 0.0) {
	       dy -= lengths_[1];
               dz -= tilt_+lengths_[2];
	       translate[1] += -1;
               translate[2] += -1;
            } else {
	       dy += lengths_[1];
               dz += tilt_+lengths_[2];
	       translate[1] += +1;
               translate[2] += +1;
            }
	    h3 = (dy*tilt_+dz*(tilt_+lengths_[2])) / d2_;
            assert(fabs(h3) <= halfd2_);
         } else {
	    done = false;
         }
      }
      }

      norm = dx*dx+dy*dy+dz*dz;
      return norm;
   }

   inline 
   double MonoclinicBoundaryMI::distanceSq(const Vector &r1, const Vector &r2) const
   {
      double dx = r1[0] - r2[0];
      double dy = r1[1] - r2[1];
      double dz = r1[2] - r2[2];
      double norm = 0.0;

      if ( fabs(dx) > halfLengths_[0] ) {
         if (dx > 0.0) {
            dx -= lengths_[0];
         } else {
            dx += lengths_[0];
         }
         assert(fabs(dx) <= halfLengths_[0]);
      }

      bool done = true;
      double h1 = 0.0;
      double h2 = 0.0;
      double h3 = 0.0; 
      
      while(done){

	 h1 = (dy*lengths_[1]+dz*tilt_) / e_;
         if ( fabs(h1) > halfe_ ) {
	    done = true;
            if (h1 > 0.0) {
               dy -= lengths_[1];
	       dz -= tilt_;
            } else {
               dy += lengths_[1];
	       dz += tilt_;
            }
            assert(fabs(dy*lengths_[1]+dz*tilt_) / e_ <= halfe_);
         } else {
	    done = false;
         }
         
         if ( fabs(dz) > halfLengths_[2] ) {
	    done = true;
            if (dz > 0.0) {
               dz -= lengths_[2];
            } else {
               dz += lengths_[2];
            }
            assert(fabs(dz) <= halfLengths_[2]);
         } else {
	    done = false;
         }

	 if ( (tilt_ > 0.0) ) {
         h2 = (dy*tilt_+dz*(tilt_-lengths_[2])) / d1_; 
         if ( fabs(h2) > halfd1_ ) {
	    done = true;
            if (h2 > 0.0) {
	       dy -= lengths_[1];
               dz -= tilt_-lengths_[2];
            } else {
	       dy += lengths_[1];
               dz += tilt_-lengths_[2];
            }
	    h2 = (dy*tilt_+dz*(tilt_-lengths_[2])) / d1_;
            assert(fabs(h2) <= halfd1_);
         } else {
	    done = false;
         }
         } else {
         h3 = (dy*tilt_+dz*(tilt_+lengths_[2])) / d2_;
         if ( fabs(h3) > halfd2_ ) {
	    done = true;
            if (h3 > 0.0) {
	       dy -= lengths_[1];
               dz -= tilt_+lengths_[2];
            } else {
	       dy += lengths_[1];
               dz += tilt_+lengths_[2];
            }
	    h3 = (dy*tilt_+dz*(tilt_+lengths_[2])) / d2_;
            assert(fabs(h3) <= halfd2_);
         } else {
	    done = false;
         }
      }
      }

      norm = dx*dx+dy*dy+dz*dz;
      return norm;
   }

   inline 
   double MonoclinicBoundaryMI::distanceSq(const Vector &r1, const Vector &r2, 
                                          Vector &dr) const
   {
      double dx = r1[0] - r2[0];
      double dy = r1[1] - r2[1];
      double dz = r1[2] - r2[2];
      double norm = 0.0;

      if ( fabs(dx) > halfLengths_[0] ) {
         if (dx > 0.0) {
            dx -= lengths_[0];
         } else {
            dx += lengths_[0];
         }
         assert(fabs(dx) <= halfLengths_[0]);
      }

      bool done = true;
      double h1 = 0.0;
      double h2 = 0.0;
      double h3 = 0.0; 
      
      while(done){

	 h1 = (dy*lengths_[1]+dz*tilt_) / e_;
         if ( fabs(h1) > halfe_ ) {
	    done = true;
            if (h1 > 0.0) {
               dy -= lengths_[1];
	       dz -= tilt_;
            } else {
               dy += lengths_[1];
	       dz += tilt_;
            }
            assert(fabs(dy*lengths_[1]+dz*tilt_) / e_ <= halfe_);
         } else {
	    done = false;
         }
         
         if ( fabs(dz) > halfLengths_[2] ) {
	    done = true;
            if (dz > 0.0) {
               dz -= lengths_[2];
            } else {
               dz += lengths_[2];
            }
            assert(fabs(dz) <= halfLengths_[2]);
         } else {
	    done = false;
         }

	 if ( (tilt_ > 0.0) ) {
         h2 = (dy*tilt_+dz*(tilt_-lengths_[2])) / d1_; 
         if ( fabs(h2) > halfd1_ ) {
	    done = true;
            if (h2 > 0.0) {
	       dy -= lengths_[1];
               dz -= tilt_-lengths_[2];
            } else {
	       dy += lengths_[1];
               dz += tilt_-lengths_[2];
            }
	    h2 = (dy*tilt_+dz*(tilt_-lengths_[2])) / d1_;
            assert(fabs(h2) <= halfd1_);
         } else {
	    done = false;
         }
         } else {
         h3 = (dy*tilt_+dz*(tilt_+lengths_[2])) / d2_;
         if ( fabs(h3) > halfd2_ ) {
	    done = true;
            if (h3 > 0.0) {
	       dy -= lengths_[1];
               dz -= tilt_+lengths_[2];
            } else {
	       dy += lengths_[1];
               dz += tilt_+lengths_[2];
            }
	    h3 = (dy*tilt_+dz*(tilt_+lengths_[2])) / d2_;
            assert(fabs(h3) <= halfd2_);
         } else {
	    done = false;
         }
      }
      }

      dr[0]= dx;
      dr[1]= dy;
      dr[2]= dz;

      norm = dx*dx+dy*dy+dz*dz;
      return norm;
   }

   /**
   * istream extractor for a MonoclinicBoundaryMI.
   *
   * \param  in        input stream
   * \param  boundary  MonoclinicBoundaryMI to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, MonoclinicBoundaryMI& boundary);

   /**
   * ostream inserter for an MonoclinicBoundaryMI.
   *
   * \param  out      output stream
   * \param  boundary MonoclinicBoundaryMI to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, 
                              const MonoclinicBoundaryMI& boundary);

   /**
   * Return actual lattice system.
   */
   inline LatticeSystem MonoclinicBoundaryMI::latticeSystem()
   { return lattice_; }

}
 
#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiTraits.h>

namespace Util
{

   /**
   * Send an MonoclinicBoundaryMI via MPI.
   */
   template <>
   void send<Util::MonoclinicBoundaryMI>(MPI::Comm& comm, 
             Util::MonoclinicBoundaryMI& data, int dest, int tag);

   /**
   * Receive an MonoclinicBoundaryMI via MPI.
   */
   template <>
   void recv<Util::MonoclinicBoundaryMI>(MPI::Comm& comm, 
             Util::MonoclinicBoundaryMI& data, int source, int tag);

   /**
   * Broadcast an MonoclinicBoundaryMI via MPI.
   */
   template <>
   void bcast<Util::MonoclinicBoundaryMI>(MPI::Intracomm& comm, 
              Util::MonoclinicBoundaryMI& data, int root);

   /**
   * Explicit specialization MpiTraits<MonoclinicBoundaryMI>.
   */
   template <>
   class MpiTraits<Util::MonoclinicBoundaryMI>
   {
   public:
      static MPI::Datatype type;         ///< MPI Datatype
      static bool hasType;               ///< Is the MPI type initialized?
   };

}
#endif // ifdef  UTIL_MPI
#endif // ifndef MONOCLINIC_BOUNDARY_H
