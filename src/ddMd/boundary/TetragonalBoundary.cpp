#ifndef TETRAGONAL_BOUNDARY_CPP
#define TETRAGONAL_BOUNDARY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TetragonalBoundary.h"
#include "LatticeSystem.h"
#include <util/space/Dimension.h>
#include <util/math/feq.h>
#include <util/format/Dbl.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /* 
   * Constructor. 
   */
   TetragonalBoundary::TetragonalBoundary() 
    : OrthoBoundaryBase()
   {}

   /* 
   * Set box length.
   */
   void TetragonalBoundary::setLengths(double ab, double c) 
   {  
      maxima_[0] = ab;
      maxima_[1] = ab;
      maxima_[2] = c;
      reset();
   }

   /* 
   * Check consistency of data.
   */
   bool TetragonalBoundary::isValid() 
   {  
      OrthoRegion::isValid(); 
      for (int i = 0; i < Dimension; ++i) {
         if (!feq(minima_[i], 0.0))
            UTIL_THROW("minima_[i] != 0");
      }
      if (!feq(lengths_[0], lengths_[1]))
         UTIL_THROW("lengths_[0] != lengths_[1]");
      return true;
   }

   /* 
   * Input a TetragonalBoundary from an istream.
   */
   std::istream& operator>>(std::istream& in, TetragonalBoundary& boundary)
   {
      LatticeSystem lattice;
      in >> lattice;
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
         UTIL_THROW("Lattice system must be tetragonal or cubic");
      }
      boundary.reset();
      return in;
   }

   /*
   * Output a TetragonalBoundary to an ostream.
   */
   std::ostream& operator<<(std::ostream& out, const TetragonalBoundary& boundary) 
   {
      LatticeSystem lattice = Tetragonal;
      out <<  lattice << "   ";
      out << Dbl(boundary.maxima_[0],20, 10);
      out << Dbl(boundary.maxima_[2],20, 10);
      return out;
   }

} 
#endif
