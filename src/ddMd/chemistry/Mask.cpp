#ifndef MASK_CPP
#define MASK_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include "Atom.h"
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   Mask::Mask()
    : size_(0)
   {
      for (int i=0; i < Capacity; ++i) {
         atomIds_[i] = -1;
      }
   }

   /*
   * Clear the mask, i.e., remove all Atoms.
   */
   void Mask::clear()
   {
      for (int i=0; i < Capacity; ++i) {
         atomIds_[i] = -1;
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
      if (isMasked(atom.id())) {
         UTIL_THROW("Attempt to add an atom to a Mask twice");
      }
      atomIds_[size_] = atom.id();
      ++size_;
   }

} 
#endif
