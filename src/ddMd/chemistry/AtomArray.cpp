#ifndef DDMD_ATOM_ARRAY_CPP
#define DDMD_ATOM_ARRAY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomArray.h"
#include "Atom.h"
#include <util/misc/Memory.h>

#include <stdlib.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   *
   * The data_ and capacity_ are nullified in Array<Atom> constructor.
   */
   AtomArray::AtomArray() 
    : Array<Atom>(),
      velocities_(0),
      masks_(0),
      plans_(0),
      ids_(0),
      contexts_(0)
   {}

   /*
   * Destructor.
   */
   AtomArray::~AtomArray()
   {
      if (data_) {
         Memory::deallocate<Atom>(data_, capacity_);
         Memory::deallocate<Vector>(velocities_, capacity_);
         Memory::deallocate<Mask>(masks_, capacity_);
         Memory::deallocate<Plan>(plans_, capacity_);
         Memory::deallocate<int>(ids_, capacity_);
         if (contexts_) {
            Memory::deallocate<AtomContext>(contexts_, capacity_);
         }
         capacity_ = 0;
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
         Log::file() << "Warning: sizeof(Atom) != 64" << std::endl;
         Log::file() << "Size of Atom  = " << sizeof(Atom)  << std::endl;
         Log::file() << "Size of Atom* = " << sizeof(Atom*) << std::endl;
      }

      // Allocate memory
      //posix_memalign((void**) &data_, 64, capacity*sizeof(Atom));
      Memory::allocate<Atom>(data_, capacity);
      Memory::allocate<Vector>(velocities_, capacity);
      Memory::allocate<Mask>(masks_, capacity);
      Memory::allocate<Plan>(plans_, capacity);
      Memory::allocate<int>(ids_, capacity);
      if (Atom::hasAtomContext()) {
         Memory::allocate<AtomContext>(contexts_, capacity);
      }
      capacity_ = capacity;

      // Initialize values.
      for (int i = 0; i < capacity_; ++i) {
        data_[i].localId_ = (i << 1);
        data_[i].arrayPtr_ = this;
        ids_[i] = -1;
        masks_[i].clear();
        plans_[i].clearFlags();
        if (Atom::hasAtomContext()) {
           contexts_[i].clear();
        }
      }

   }

   /*
   * Set all forces to zero.
   */
   void AtomArray::zeroForces() 
   {
      for (int i = 0; i < capacity_; ++i) {
        data_[i].force_[0] = 0.0;
        data_[i].force_[1] = 0.0;
        data_[i].force_[2] = 0.0;
      }
   }

   /*
   * Return true if this is already allocated, false otherwise.
   */
   bool AtomArray::isAllocated() const 
   {  return (bool)data_; }

}
#endif
