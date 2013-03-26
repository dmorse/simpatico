#ifndef DDMD_ATOM_ARRAY_CPP
#define DDMD_ATOM_ARRAY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomArray.h"
#include "Atom.h"

#include <stdlib.h>

namespace DdMd
{

   /*
   * Constructor.
   *
   * The data_ and capacity_ are nullified in Array<Atom> constructor.
   */
   AtomArray::AtomArray() 
    : Array<Atom>(),
      velocities_(0),
      plans_(0),
      masks_(0),
      ids_(0)
   {}

   /*
   * Destructor.
   */
   AtomArray::~AtomArray()
   {
      if (data_) {
         delete [] data_;
         //free(data_);
         delete [] velocities_;
         delete [] masks_;
         delete [] plans_;
         delete [] ids_;
      }
   }

   /*
   * Allocate the underlying C array.
   */
   void AtomArray::allocate(int capacity) 
   {
      if (!(data_ == 0)) {
         UTIL_THROW("Cannot re-allocate an AtomArray");
      }
      if (capacity <= 0) {
         UTIL_THROW("Cannot allocate with capacity <= 0");
      }
      if (sizeof(Atom) != 64) {
         std::cout << "Warning: sizeof(Atom) != 64" << std::endl;
         std::cout << "Size of Atom  = " << sizeof(Atom)  << std::endl;
         std::cout << "Size of Atom* = " << sizeof(Atom*) << std::endl;
      }

      // Allocate memory
      //posix_memalign((void**) &data_, 64, capacity*sizeof(Atom));
      data_        = new Atom[capacity];
      velocities_  = new Vector[capacity];
      masks_       = new Mask[capacity];
      plans_       = new Plan[capacity];
      ids_         = new int[capacity];
      capacity_    = capacity;

      // Initialize values.
      unsigned int localId;
      for (int i=0; i < capacity; ++i) {
        data_[i].localId_ = (i << 1);
        data_[i].arrayPtr_ = this;
        localId = (data_[i].localId_ >> 1);
        ids_[i] = -1;
        masks_[i].clear();
        plans_[i].clearFlags();
      }

   }

   /*
   * Return true if this is already allocated, false otherwise.
   */
   bool AtomArray::isAllocated() const 
   {  return !(data_ == 0); }

}
#endif
