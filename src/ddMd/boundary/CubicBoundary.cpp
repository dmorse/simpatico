#ifndef CUBIC_BOUNDARY_CPP
#define CUBIC_BOUNDARY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CubicBoundary.h"
#include "LatticeSystem.h"
#include <util/format/Dbl.h>
#include <util/space/Dimension.h>
#include <util/math/feq.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /* 
   * Constructor. 
   */
   CubicBoundary::CubicBoundary() 
    : OrthoBoundaryBase()
   {}

   /* 
   * Set box length.
   */
   void CubicBoundary::setLength(double length) 
   {  
      assert(length > 0.0);
      maxima_[0] = length;
      maxima_[1] = length;
      maxima_[2] = length;
      reset();
   }

   /* 
   * Return true if internal state is valid, or throw Exception.
   */
   bool CubicBoundary::isValid() 
   {  
      OrthoRegion::isValid(); 
      for (int i = 0; i < Dimension; ++i) {
         if (!feq(minima_[i], 0.0))
            UTIL_THROW("minima_[i] != 0");
      }
      if (!feq(lengths_[0], lengths_[1]))
         UTIL_THROW("lengths_[0] != lengths_[1]");
      if (!feq(lengths_[1], lengths_[2]))
         UTIL_THROW("lengths_[1] != lengths_[2]");
      return true;
   }

   /* 
   * Input a CubicBoundary from an istream.
   */
   std::istream& operator>>(std::istream& in, CubicBoundary& boundary)
   {
      LatticeSystem lattice;
      double        length;

      in >> lattice;
      if (lattice != Cubic) 
         UTIL_THROW("Lattice type must be Cubic");

      in >> length;
      if (length <= 0.0) 
         UTIL_THROW("Input length <= 0.0");

      boundary.maxima_[0] = length;
      boundary.maxima_[1] = length;
      boundary.maxima_[2] = length;
      boundary.reset();
      return in;
   }

   /* 
   * Output a CubicBoundary to an ostream.
   */
   std::ostream& operator<<(std::ostream& out, const CubicBoundary& boundary) 
   {
      LatticeSystem lattice = Cubic;
      out << lattice << "    " << Dbl(boundary.maxima_[0], 20, 10);
      return out;
   }

} 
#endif
