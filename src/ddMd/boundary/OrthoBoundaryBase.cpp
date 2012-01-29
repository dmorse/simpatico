#ifndef ORTHO_BOUNDARY_BASE_CPP
#define ORTHO_BOUNDARY_BASE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthoBoundaryBase.h"
#include <util/random/Random.h>
#include <util/math/Constants.h>

namespace DdMd
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   OrthoBoundaryBase::OrthoBoundaryBase() 
    : OrthoRegion()
   {
      for (int i = 0; i < Dimension; ++i) {

         bravaisBasisVectors_.append(Vector::Zero);
         bravaisBasisVectors_[i][i] = lengths_[i];

         reciprocalBasisVectors_.append(Vector::Zero);
         reciprocalBasisVectors_[i][i] = 2.0*Constants::Pi/lengths_[i];

      }
   }

   /* 
   * Generate a random position within the box
   *
   * \param random random number generator object
   * \param r[3]   vector of random coordinates
   */
   void OrthoBoundaryBase::randomPosition(Random &random, Vector &r) const 
   {
     for (int i=0; i < Dimension; ++i) {
        r[i] = random.getFloat(minima_[i], maxima_[i]);
     }
   }

   /* 
   * Reset all quantities that depend on unit cell lengths.
   */
   void OrthoBoundaryBase::reset()
   {
      resetRegion();
      for (int i = 0; i < Dimension; ++i) {
         bravaisBasisVectors_[i][i] = lengths_[i];
         reciprocalBasisVectors_[i][i] = 2.0*Constants::Pi/lengths_[i];
      }
   }

} 
#endif
