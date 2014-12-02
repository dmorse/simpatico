/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   Mask::Mask()
    : size_(0)
   {
      for (int i=0; i < Capacity; ++i) {
         atomPtrs_[i] = 0;
      }
   }

   /*
   * Clear the mask, i.e., remove all Atoms.
   */
   void Mask::clear()
   {
      for (int i=0; i < Capacity; ++i) {
         atomPtrs_[i] = 0;
      }
      size_ = 0;
   }

   /*
   * Add an an atom to the mask.
   */
   void Mask::append(const Atom& atom)
   {
      if (size_ >= Capacity) {
         UTIL_THROW("Too many masked partners for one Atom");
      }
      if (isMasked(atom)) {
         UTIL_THROW("Attempt to add an atom to a Mask twice");
      }
      atomPtrs_[size_] = &atom;
      ++size_;
   }

} 
