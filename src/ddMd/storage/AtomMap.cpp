#ifndef DDMD_ATOM_MAP_CPP
#define DDMD_ATOM_MAP_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomMap.h"
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   AtomMap::AtomMap()
    : atomPtrs_(),
      nLocal_(0),
      nGhost_(0),
      totalAtomCapacity_(0),
      isInitialized_(false)
   {}
 
   /*
   * Destructor.
   */
   AtomMap::~AtomMap()
   {}

   /*
   * Allocate and initialize all containers (private).
   */
   void AtomMap::allocate(int totalAtomCapacity)
   {
      // Precondition
      if (isInitialized_) {
         UTIL_THROW("AtomMap can only be initialized once");
      }
      atomPtrs_.allocate(totalAtomCapacity);
      totalAtomCapacity_ = totalAtomCapacity;
      for (int i = 0; i < totalAtomCapacity_; ++i) {
         atomPtrs_[i] = 0;
      }
      isInitialized_ = true;
   }

   /*
   * Register new local Atom in internal data structures.
   */ 
   void AtomMap::addLocal(Atom* ptr)
   {
      int atomId = ptr->id();

      // Preconditions
      if (atomId < 0 || atomId >= totalAtomCapacity_) {
         std::cout << "atomId = " << atomId << std::endl;
         UTIL_THROW("atomId is out of range");
      }
      if (find(atomId)) {
         std::cout << "atomId       = " << atomId << std::endl;
         std::cout << "New Position = " << ptr->position() 
                   << std::endl;
         std::cout << "Old Position = " << find(atomId)->position() 
                   << std::endl;
         UTIL_THROW("Atom with specified id is already present");
      }

      // Add 
      atomPtrs_[atomId] = ptr;
      ++nLocal_;
   }

   /*
   * Remove a specific local Atom.
   */
   void AtomMap::removeLocal(Atom* ptr)
   {
      int atomId = ptr->id();
      if (atomId < 0 || atomId >= totalAtomCapacity_) {
         std::cout << "atomId = " << atomId << std::endl;
         UTIL_THROW("atomId is out of range");
      }
      if (ptr == atomPtrs_[atomId]) {
         atomPtrs_[atomId] = 0;
         --nLocal_;
      } else {
         if (0 == atomPtrs_[atomId]) {
            UTIL_THROW("Error: Attempt to remove absent local atom");
         } else {
            UTIL_THROW("Error: Inconsistent pointer");
         }
      }
   }

   // Ghost atom mutators

   /*
   * Register new ghost Atom to internal data structures.
   */ 
   void AtomMap::addGhost(Atom* ptr)
   {
      int atomId = ptr->id();
      if (atomId < 0 || atomId >= totalAtomCapacity_) {
         UTIL_THROW("atomId is out of range");
      }

      #if 0
      // Add iff no particle with same id is already present
      if (0 == atomPtrs_[atomId]) {
         atomPtrs_[atomId] = ptr;
         ++nGhost_;
      }
      #endif

      std::pair<GhostMap::iterator, bool> ret;
      ret = ghostMap_.insert(std::pair<int, Atom*>(atomId, ptr));
      if (ret.second) {
         ++nGhost_;
      }

   }

   /*
   * Remove a specific ghost Atom.
   */
   void AtomMap::removeGhost(Atom* ptr)
   {
      int atomId = ptr->id();
      if (atomId < 0 || atomId >= totalAtomCapacity_) {
         UTIL_THROW("atomId is out of range");
      }

      #if 0
      if (atomPtrs_[atomId] == ptr) {
         atomPtrs_[atomId] = 0;
         --nGhost_;
      } else {
         if (0 == atomPtrs_[atomId]) {
            UTIL_THROW("Error: Attempt to remove absent atom");
         }
         // Note: No error if another atom with same id is present.
      }
      #endif

      GhostMap::iterator iter;
      iter = ghostMap_.find(atomId);
      if (iter != ghostMap_.end()) {
         if (iter->second == ptr) {
            ghostMap_.erase(iter);
            --nGhost_;
         }
         // Note: No error if another atom with same id is present.
      } else {
         UTIL_THROW("Error: Attempt to remove absent ghost");
      }

   }

   /*
   * Check validity of this AtomMap.
   *
   * Returns true if all is ok, or throws an Exception.
   */
   bool AtomMap::isValid() const
   {
      Atom* ptr;
      int i, j;
      j = 0;
      for (i = 0; i < totalAtomCapacity_ ; ++i) {
         ptr = atomPtrs_[i];
         if (ptr != 0) {
            ++j;
            if (ptr->id() != i) {
               std::cout << std::endl;
               std::cout << "Index i in atomPtrs_  " << i << std::endl;
               std::cout << "atomPtrs_[i]->id()    " << ptr->id() << std::endl;
               UTIL_THROW("ptr->id() != i");
            }
         }
      }
      if (j != nLocal_) {
         UTIL_THROW("Inconsistent count of local atoms in AtomMap");
      }
      if (ghostMap_.size() != nGhost_) {
         UTIL_THROW("Inconsistent count of ghost atoms in AtomMap");
      }
      return true;
   }

}
#endif
