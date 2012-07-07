#ifndef MONOCLINIC_BOUNDARY_CPP
#define MONOCLINIC_BOUNDARY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MonoclinicBoundary.h"
#include <util/space/Vector.h>
#include <util/space/Dimension.h>
#include <util/random/Random.h>
#include <util/math/Constants.h>
#include <util/math/feq.h>
#include <util/format/Dbl.h>
#include <util/global.h>

namespace Util
{

   /* 
   * Default constructor.
   */
   MonoclinicBoundary::MonoclinicBoundary() 
    : lattice_(Monoclinic)
   {
      for (int i = 0; i < Dimension; ++i) {

	 MonoclinicBoundary::minima_[i] = 0.0;
	 MonoclinicBoundary::maxima_[i] = 1.0;
	 MonoclinicBoundary::lengths_[i] = 1.0;
	 MonoclinicBoundary::halfLengths_[i] = 0.5;

         bravaisBasisVectors_.append(Vector::Zero);
         bravaisBasisVectors_[i][i] = lengths_[i];
         reciprocalBasisVectors_.append(Vector::Zero);
         reciprocalBasisVectors_[i][i] = 2.0 * Constants::Pi / lengths_[i];

      }

	 bravaisBasisVectors_[1][2] = MonoclinicBoundary::tilt_;
         reciprocalBasisVectors_[2][1] = -2.0 * Constants::Pi*(MonoclinicBoundary::tilt_/(lengths_[1]*lengths_[2]));

	 MonoclinicBoundary::tilt_ = 0.0;
	 MonoclinicBoundary::volume_ = 1.0;
	 
	 MonoclinicBoundary::c1_ = 1.0;
	 MonoclinicBoundary::c2_ = 1.0;
	 MonoclinicBoundary::c3_ = 0.0;
	 MonoclinicBoundary::e_ = 1.0;
	 MonoclinicBoundary::halfe_ = 0.5;
	 MonoclinicBoundary::maxValidityDistanceSq_ = 0.5;
   }

   /* 
   * Set box lengths and then call reset.
   */
   void MonoclinicBoundary::setLengths(const Vector &lengths, const double d) 
   {  
      maxima_  = lengths; 
      tilt_ = d;
      lattice_ = Monoclinic;
      reset(); 
   }

   /* 
   * Reset all quantities that depend on unit cell lengths.
   */
   void MonoclinicBoundary::reset()
   {
      for (int i = 0; i < Dimension; ++i) {
	 assert(maxima_[i] > minima_[i]);
         lengths_[i] = maxima_[i] - minima_[i];
         halfLengths_[i] = 0.5*lengths_[i];
         bravaisBasisVectors_[i][i] = lengths_[i];
         reciprocalBasisVectors_[i][i] = 2.0*Constants::Pi/lengths_[i];
      }
	 bravaisBasisVectors_[1][2] = MonoclinicBoundary::tilt_;
         reciprocalBasisVectors_[2][1] = -2.0*Constants::Pi*(MonoclinicBoundary::tilt_/(lengths_[1]*lengths_[2]));

         volume_ = lengths_[0] * lengths_[1] * lengths_[2];

	 e_ = sqrt(lengths_[1]*lengths_[1]+tilt_*tilt_); 	 
	 halfe_ = sqrt(lengths_[1]*lengths_[1]+tilt_*tilt_) / 2.0; 	 
	 c1_ = sqrt(lengths_[1]*lengths_[1]+tilt_*tilt_) / lengths_[1];
	 c2_ = 1.0;
	 c3_ = -tilt_/(lengths_[1]);
	 maxValidityDistanceSq_ = maxValidityDistanceSq();
   }

   /* 
   * Generate a random position within the box and then tilt_ the box.
   *
   * \param random random number generator object
   * \param r      Vector of random coordinates
   */
   void MonoclinicBoundary::randomPosition(Random &random, Vector &r) const 
   {
     for (int i=0; i < Dimension; ++i) {
        r[i] = random.uniform(minima_[i], maxima_[i]);
     }
   }

   /* 
   * Check consistency of data.
   */
   bool MonoclinicBoundary::isValid() 
   {  
      for (int i = 0; i < Dimension; ++i) {
         if (maxima_[i] <= minima_[i])   
            UTIL_THROW("maxima_[i] <= minima_[i]");
         if (!feq(lengths_[i], maxima_[i] - minima_[i]))
            UTIL_THROW("lengths_[i] != maxima_[i] - minima_[i]");
         if (!feq(halfLengths_[i], 0.5*lengths_[i]))
            UTIL_THROW("halfLengths_[i] != 0.5*lengths_[i]");
         if (!feq(minima_[i], 0.0))
            UTIL_THROW("minima_[i] != 0");
      }
      if (!feq(volume_, lengths_[0]*lengths_[1]*lengths_[2]))
         UTIL_THROW("volume_ != product of lengths_");
      return true;
   }

   /* 
   * Input a MonoclinicBoundary from an istream, without line breaks.
   */
   std::istream& operator>>(std::istream& in, MonoclinicBoundary &boundary)
   {
      LatticeSystem lattice;
      in >> lattice;
      if (lattice == Monoclinic) {
         in >> boundary.maxima_;
      } else {
         UTIL_THROW("Lattice must be Monoclinic");
      }
      boundary.lattice_ = lattice;
      boundary.reset();
      return in;
   }

   /* 
   * Output an MonoclinicBoundary to an ostream, without line breaks.
   */
   std::ostream& 
   operator<<(std::ostream& out, const MonoclinicBoundary &boundary) 
   {
      out << boundary.lattice_ << "   ";
      if (boundary.lattice_ == Monoclinic) {
         out << boundary.lengths_;
      } 

      return out;
   }
 
   #ifdef UTIL_MPI
   template <>
   void send<Util::MonoclinicBoundary>(MPI::Comm& comm, 
             Util::MonoclinicBoundary& data, int dest, int tag)
   {
      Vector lengths = data.lengths();
      send<Vector>(comm, lengths, dest, tag);
   }

   template <>
   void recv<Util::MonoclinicBoundary>(MPI::Comm& comm, 
             Util::MonoclinicBoundary& data, int source, int tag)
   {
      Vector lengths;
      recv<Vector>(comm, lengths, source, tag);
      data.setLengths(lengths);
   }

   template <>
   void bcast<Util::MonoclinicBoundary>(MPI::Intracomm& comm, 
              Util::MonoclinicBoundary& data, int root)
   {
      Vector lengths; 
      int    rank = comm.Get_rank();
      if (rank == root) 
         lengths = data.lengths();
      bcast<Vector>(comm, lengths, root);
      if (rank != root) 
         data.setLengths(lengths);
   }

   /*
   * Initialize MPI Datatype.
   */
   MPI::Datatype MpiTraits<Util::MonoclinicBoundary>::type = MPI::BYTE;
   bool MpiTraits<Util::MonoclinicBoundary>::hasType = false;
   #endif

}

#endif
