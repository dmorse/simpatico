#ifndef SPAN_ATOM_STORAGE_CPP
#define SPAN_ATOM_STORAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomStorage.h"

namespace SpAn 
{

   /*
   * Constructor.
   */
   AtomStorage::AtomStorage()
    : newPtr_(0)
   {}

   /*
   * Destructor.
   */
   AtomStorage::~AtomStorage()
   {}

   /*
   * Allocate and initialize memory.
   */
   void AtomStorage::allocate(int capacity)
   {
      atoms_.allocate(capacity);
      atomPtrs_.allocate(capacity);
      clear();
   }

   /*
   * Return pointer to location for new atom.
   */
   Atom* AtomStorage::newPtr()
   {
      if (newPtr_) {
         UTIL_THROW("Error: an new atom is still active");
      }
      int size = atoms_.size() + 1;
      atoms_.resize(size);
      newPtr_ = &atoms_[size - 1];
      return newPtr_;
   }

   /*
   * Finalize addition of new atom.
   */
   void AtomStorage::add()
   {
      if (!newPtr_) {
         UTIL_THROW("Error: No active new atom");
      }
      int id = newPtr_->id;
      atomPtrs_[id] = newPtr_;
      newPtr_ = 0;
   }

   /*
   * Remove all atoms and bonds - set to empty state.
   */
   void AtomStorage::clear()
   {
      atoms_.clear();
      int capacity = atomPtrs_.capacity();
      for (int i = 0; i < capacity; ++i) {
         atomPtrs_[i] = 0;
      }
   }

}
#endif
