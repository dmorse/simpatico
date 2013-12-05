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

//#define DDMD_GROUP_STORAGE_DEBUG

namespace DdMd
{
 
   using namespace Util;

   /*
   * Find pointers to atoms in Group<2> after exchanging all ghosts.
   */
   template <>
   void GroupStorage<2>::findGhosts(AtomStorage& atomStorage, 
                                    const Boundary& boundary)
   {
      const AtomMap& atomMap = atomStorage.map();
      GroupIterator<2> groupIter;
      Atom* aPtr;  // Pointer to atom a (must be local)
      Atom* bPtr;  // Pointer to atom b (may be local or ghost)
      int k;       // id for atom b within group
      for (begin(groupIter); groupIter.notEnd(); ++groupIter) {
         aPtr = groupIter->atomPtr(0);
         if (aPtr) {
            k = 1;
         } else {
            aPtr = groupIter->atomPtr(1);
            assert(aPtr);
            k = 0;
         }
         atomMap.findNearestImage(groupIter->atomId(k), aPtr->position(), 
                                  boundary, bPtr);
         groupIter->setAtomPtr(k, bPtr);
      }
   }

} // namespace DdMd
#endif
