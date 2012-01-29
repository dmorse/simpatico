#ifndef ORTHORHOMBIC_BOUNDARY_CPP
#define ORTHORHOMBIC_BOUNDARY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthorhombicBoundary.h"
#include "LatticeSystem.h"
#include <util/space/Vector.h>
#include <util/space/Dimension.h>
#include <util/math/feq.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /* 
   * Constructor. Zeros all data members.
   */
   OrthorhombicBoundary::OrthorhombicBoundary() 
    : OrthoBoundaryBase()
   {}

   /* 
   * Set box length array lengths_[], and then call initialize().
   */
   void OrthorhombicBoundary::setLengths(const Vector &lengths) 
   {  
      maxima_ = lengths; 
      reset(); 
   }

   /* 
   * Check consistency of data.
   */
   bool OrthorhombicBoundary::isValid() 
   {  
      OrthoRegion::isValid(); 
      for (int i = 0; i < Dimension; ++i) {
         if (!feq(minima_[i], 0.0))
            UTIL_THROW("minima_[i] != 0");
      }
      return true;
   }

   /* 
   * Input a OrthorhombicBoundary from an istream, without line breaks.
   */
   std::istream& operator>>(std::istream& in, OrthorhombicBoundary &boundary)
   {
      LatticeSystem lattice;
      in >> lattice;
      if (lattice == Orthorhombic) {
         in >> boundary.maxima_;
      } else
      if (lattice == Tetragonal) {
         double ab, c;
         in >> ab >> c;
         boundary.maxima_[0] = ab;
         boundary.maxima_[1] = ab;
         boundary.maxima_[2] = c;
      } else 
      if (lattice == Cubic) {
         double a; 
         in >> a; 
         boundary.maxima_[0] = a;
         boundary.maxima_[1] = a;
         boundary.maxima_[2] = a;
      } else {
         UTIL_THROW("Lattice must be orthorhombic, tetragonal or cubic");
      }
      boundary.reset();
      return in;
   }

   /* 
   * Output a OrthorhombicBoundary to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const OrthorhombicBoundary &boundary) 
   {
      LatticeSystem lattice = Orthorhombic;
      out << lattice << "   ";
      out << boundary.lengths_;
      return out;
   }

}
 
#ifdef UTIL_MPI
namespace Util
{

   template <>
   void send<DdMd::OrthorhombicBoundary>(MPI::Comm& comm, 
             DdMd::OrthorhombicBoundary& data, int dest, int tag)
   {
      Vector lengths = data.lengths();
      send<Vector>(comm, lengths, dest, tag);
   }

   template <>
   void recv<DdMd::OrthorhombicBoundary>(MPI::Comm& comm, 
             DdMd::OrthorhombicBoundary& data, int source, int tag)
   {
      Vector lengths;
      recv<Vector>(comm, lengths, source, tag);
      data.setLengths(lengths);
   }

   template <>
   void bcast<DdMd::OrthorhombicBoundary>(MPI::Intracomm& comm, 
              DdMd::OrthorhombicBoundary& data, int root)
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
   MPI::Datatype MpiTraits<DdMd::OrthorhombicBoundary>::type = MPI::BYTE;
   bool MpiTraits<DdMd::OrthorhombicBoundary>::hasType = false;

}
#endif

#endif
