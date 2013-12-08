#ifndef DDMD_GROUP_STORAGE_CPP
#define DDMD_GROUP_STORAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "GroupStorage.h"
#include "AtomStorage.h"
#include "AtomMap.h"

//#define DDMD_GROUP_STORAGE_DEBUG

namespace DdMd
{
 
   using namespace Util;

   /*
   * Find pointers to atoms in Group<2> after exchanging all ghosts.
   */
   template <>
   void GroupStorage<2>::finishGhostExchange(const AtomMap& map, 
                                             const Boundary& boundary)
   {
      GroupIterator<2> iter;
      Atom* aPtr;  // Pointer to atom a (must be local)
      Atom* bPtr;  // Pointer to atom b (may be local or ghost)
      Atom* cPtr;  // Pointer to unused local image of b atom
      int k;       // id for atom b within group
      for (begin(iter); iter.notEnd(); ++iter) {
         aPtr = iter->atomPtr(0);
         if (aPtr) {
            k = 1;
         } else {
            aPtr = iter->atomPtr(1);
            assert(aPtr);
            k = 0;
         }
         assert(!aPtr->isGhost());
         cPtr = map.findNearestImage(iter->atomId(k), aPtr->position(), 
                                     boundary, bPtr);
         if (cPtr) {
            makeGroupImage(*iter, k, cPtr, map, boundary);
         }
         iter->setAtomPtr(k, bPtr);
      }
   }

   /*
   * Make a new image of a group.
   */
   template <>
   void GroupStorage<2>::makeGroupImage(Group<2>& group, 
                                        int rootId, Atom* rootPtr,
                                        const AtomMap& map, 
                                        const Boundary& boundary)
   {
      // Get a new group object and add to group and ghost sets
      Group<2>* newPtr = &reservoir_.pop();
      groupSet_.append(*newPtr);
      ghosts_.append(*newPtr);
      if (groupSet_.size() > maxNGroupLocal_) {
         maxNGroupLocal_ = groupSet_.size();
      }

      // Copy group id, atomIds, & root pointer to new group
      newPtr->setId(group.id());
      for (int j = 0; j < 2; ++j) {
         newPtr->setAtomId(j, group.atomId(j));
         newPtr->clearAtomPtr(j);
      }
      newPtr->setAtomPtr(rootId, rootPtr);

      Atom* bPtr; // pointer to non-root atom
      int k;      // index of non-root atom within the group
      k = (rootId == 0) ? 1 : 0;
      map.findNearestImage(newPtr->atomId(k), rootPtr->position(), 
                           boundary, bPtr);
      newPtr->setAtomPtr(k, bPtr);

      // Note: Their can be atom most one image of a Group<2>, 
      // so makeGroupImage does not need to call itself for N=2.
   }

} // namespace DdMd
#endif
