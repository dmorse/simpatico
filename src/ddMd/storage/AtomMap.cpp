#ifndef DDMD_ATOM_MAP_CPP
#define DDMD_ATOM_MAP_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomMap.h"
#include <ddMd/storage/ConstGhostIterator.h>
#include <util/containers/ArraySet.h>
#include <util/global.h>

#include <utility>

namespace DdMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   AtomMap::AtomMap()
    : atomPtrs_(),
      nLocal_(0),
      nGhostDistinct_(0),
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

      // Preconditions
      int atomId = ptr->id();
      if (atomId < 0 || atomId >= totalAtomCapacity_) {
         Log::file() << "atomId = " << atomId << std::endl;
         UTIL_THROW("atomId is out of range");
      }
      if (atomPtrs_[atomId]) {
         Log::file() << "atomId       = " << atomId << std::endl;
         Log::file() << "New Position = " << ptr->position() 
                   << std::endl;
         Log::file() << "Old Position = " << atomPtrs_[atomId]->position() 
                   << std::endl;
         UTIL_THROW("Local atom with specified id is already present");
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
         Log::file() << "atomId = " << atomId << std::endl;
         UTIL_THROW("atomId is out of range");
      }
      if (ptr == atomPtrs_[atomId]) {

         // Remove from atomPtrs_
         atomPtrs_[atomId] = 0;
         --nLocal_;

         // If possible, move an atom from ghostMap to atomPtrs_
         if (!ghostMap_.empty()) {
            GhostMap::iterator iter = ghostMap_.find(atomId);
            if (iter != ghostMap_.end()) {
               atomPtrs_[atomId] = iter->second;
               ++nGhostDistinct_;
               ghostMap_.erase(iter);
            }
         }

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

      if (0 == atomPtrs_[atomId]) {
         atomPtrs_[atomId] = ptr;
         ++nGhostDistinct_;
      } else {
         ghostMap_.insert(std::pair<int, Atom*>(atomId, ptr));
         //ghostMap_.emplace(atomId, ptr);
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

      if (atomPtrs_[atomId] == ptr) {

         // Remove from atomPtrs array
         atomPtrs_[atomId] = 0;
         --nGhostDistinct_;

         // If possible, move an atom from ghostMap to atomPtrs_
         if (!ghostMap_.empty()) {
            GhostMap::iterator iter = ghostMap_.find(atomId);
            if (iter != ghostMap_.end()) {
               atomPtrs_[atomId] = iter->second;
               ++nGhostDistinct_;
               ghostMap_.erase(iter);
            }
         }

      } else { // If ptr is not found in atomPtrs

         if (atomPtrs_[atomId] != 0) { 
            // Search ghost map
            std::pair<GhostMap::iterator, GhostMap::iterator> ret;
            ret = ghostMap_.equal_range(atomId);
            GhostMap::iterator it = ret.first;
            GhostMap::iterator last = ret.second;
            for ( ; it != last; ++it) {
               assert(it->first == atomId);
               if (it->second == ptr) {
                  ghostMap_.erase(it);
                  return;
               }
            } 
            UTIL_THROW("Error: Attempt to remove absent ghost");
         } else {
            UTIL_THROW("Error: Attempt to remove absent ghost");
         }
      }

   }

   void AtomMap::clearGhosts(const ArraySet<Atom>& ghostSet)
   {
      // Precondition
      if (ghostSet.size() != nGhost()) {
         UTIL_THROW("Inconsistent ghost set sizes"); 
      }

      // Clear extra ghost images from ghostMap_
      ghostMap_.clear();

      // Clear ghosts from atomPtrs_ array
      ConstGhostIterator iter;
      const Atom* ptr;
      int id;
      for (ghostSet.begin(iter); iter.notEnd(); ++iter) {
         id = iter->id();
         ptr = iter.get();
         if (atomPtrs_[id] == ptr) {
            atomPtrs_[id] = 0;
            --nGhostDistinct_;
         }
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
      int i, id, nAtom;

      // Validate atomPtrs_
      nAtom = 0;
      for (i = 0; i < totalAtomCapacity_ ; ++i) {
         ptr = atomPtrs_[i];
         if (ptr != 0) {
            id = ptr->id();
            if (id != i) {
               Log::file() << std::endl;
               Log::file() << "Index i in atomPtrs_  " << i << std::endl;
               Log::file() << "atomPtrs_[i]->id()    " << id << std::endl;
               UTIL_THROW("ptr->id() != i");
            }
            ++nAtom;
         }
      }
      if (nAtom != nLocal_ + nGhostDistinct_) {
         UTIL_THROW("Inconsistent count of atoms in atomPtrs");
      }

      // Validate ghostMap_
      GhostMap::const_iterator it;  
      for (it = ghostMap_.begin(); it != ghostMap_.end(); ++it) {
         id  = it->first;
         ptr = it->second;
         if (id != ptr->id()) {
            Log::file() << std::endl;
            Log::file() << "key ghostMap " << id << std::endl;
            Log::file() << "Atom::id()   " << ptr->id() << std::endl;
            UTIL_THROW("Inconsistent key in ghostMap");
         }
         if (atomPtrs_[id] == 0) {
            UTIL_THROW("Id in ghostMap_ does not appear in atomPtrs_");
         } 
      }
      return true;
   }

}
#endif
