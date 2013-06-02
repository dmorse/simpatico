#ifndef DDMD_GROUP_STORAGE_TPP
#define DDMD_GROUP_STORAGE_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "GroupStorage.h"

namespace DdMd
{

   /*
   * Identify groups that span boundaries.
   *
   * This method is called by exchangeAtoms, after computing plans for
   * exchanging atoms, based on their position, but before exchanging
   * atoms, and before clearing ghosts from any previous exchange.
   *
   * Algorithm: Loop over all Group<N> objects in the storage, identify
   * groups that span boundaries of the processor domain associated with
   * each of 6 transfer directions (3 Cartesian directions, and transfer
   * "up" and "down" in each direction). This requires information about 
   * positions of ghost as well as local atoms. For each boundary of the 
   * domain, identify atoms whose positions are "inside" and "outside".
   * Count ghost atoms very near the boundary as both inside and outside,
   * for saftey. If a group has atoms both inside and outside a domain 
   * boundary, it is marked for sending in the associated communication 
   * step.
   *
   * After calculating a ghost communication plan for each group, clear 
   * the pointers to all ghost atoms in the group. The exchangeAtoms 
   * method will clear the actual ghost atoms from the AtomStorage.
   */
   template <int N> void 
   GroupStorage<N>::markSpanningGroups(FMatrix<double, Dimension, 2>& bound, 
                                       FMatrix<double, Dimension, 2>& inner, 
                                       FMatrix<double, Dimension, 2>& outer, 
                                       IntVector& gridFlags)
   {
      double coordinate;
      GroupIterator<N> groupIter;
      Atom* atomPtr;
      int nIn;
      int nOut;
      int i, j, k;
      bool choose;

      // Loop over groups
      begin(groupIter);
      for ( ; groupIter.notEnd(); ++groupIter) {

         // Compute ghost communication plan for group
         groupIter->plan().clearFlags();
         for (i = 0; i < Dimension; ++i) {
            if (gridFlags[i]) {
               for (j = 0; j < 2; ++j) {
                  choose = false;
                  nIn = 0;
                  nOut = 0;
                  // Loop over atoms in group
                  for (k = 0; k < N; ++k) {
                     atomPtr = groupIter->atomPtr(k);
                     if (atomPtr) {
                        coordinate = atomPtr->position()[i];
                        if (atomPtr->isGhost()) {
                           if (j == 0) {
                              assert(inner(i, j) > bound(i, j));
                              if (coordinate < inner(i, j)) {
                                 ++nOut;
                              }
                              if (coordinate > outer(i, j)) {
                                 ++nIn;
                              }
                           } else { // if j = 1
                              assert(inner(i, j) < bound(i, j));
                              if (coordinate > inner(i, j)) {
                                 ++nOut;
                              }
                              if (coordinate < outer(i, j)) {
                                 ++nIn;
                              }
                           }
                        } else { // if atomPtr points to local atom
                           if (atomPtr->plan().exchange(i, j)) {
                              ++nOut;
                           } else {
                              ++nIn;
                           }
                        }
                     } else { // if atomPtr is null
                        choose = true;
                        break;
                     }
                  } // end for k (atoms in group)
                  if (nOut > 0 && nIn > 0) {
                     choose = true;
                  }
                  if (choose) {
                     groupIter->plan().setGhost(i, j);
                  } else {
                     groupIter->plan().clearGhost(i, j);
                  }
               } // end for j = 0, 1

            } // end if gridFlags[i]
         } // end for i (Cartesian axes)

         // Clear pointers to all ghost atoms in this group
         for (k = 0; k < N; ++k) {
            atomPtr = groupIter->atomPtr(k);
            if (atomPtr) {
               if (atomPtr->isGhost()) {
                  groupIter->clearAtomPtr(k);
               }
            }
         }

      }
   }

   #ifdef UTIL_MPI
   /*
   * Pack groups for exchange.
   *
   * Pack groups that contain atoms marked for exchange in this 
   * direction (direction i, j), and remove empty groups.
   *
   * Algorithm: Loop over groups. If the group contains one or
   * more atoms that are marked for exchange in direction i, j, 
   * pack the group for sending along. Remove empty groups in
   * a separate loop.
   */
   template <int N>
   void GroupStorage<N>::pack(int i, int j, Buffer& buffer)
   {
      GroupIterator<N> groupIter;
      Atom* atomPtr;
      int k, nAtom;
      bool choose;
      emptyGroups_.clear();

      // Pack Groups
      buffer.beginSendBlock(Buffer::GROUP2 + N - 2);
      begin(groupIter);
      for ( ; groupIter.notEnd(); ++groupIter) {
         choose = false;
         nAtom = 0;
         for (k = 0; k < N; ++k) {
            atomPtr = groupIter->atomPtr(k);
            if (atomPtr) {
               if (atomPtr->plan().exchange(i, j)) {
                  choose = true;
                  groupIter->clearAtomPtr(k);
               } else {
                  ++nAtom;
               }
            }
         }
         if (nAtom == 0) {
            emptyGroups_.append(*groupIter);
         }
         if (choose) {
            groupIter->pack(buffer);
         }
      }
      buffer.endSendBlock();

      // Remove empty groups
      int nEmpty = emptyGroups_.size();
      for (int k = 0; k < nEmpty; ++k) {
         remove(&(emptyGroups_[k]));
      }
   }

   /*
   * Unpack groups into GroupStorage.
   */
   template <int N>
   void GroupStorage<N>::unpack(Buffer& buffer, AtomStorage& atomStorage)
   {
      Group<N>* newGroupPtr;
      Group<N>* oldGroupPtr;
      int groupId;

      buffer.beginRecvBlock();
      while (buffer.recvSize() > 0) {
         newGroupPtr = newPtr();
         newGroupPtr->unpack(buffer);
         groupId = newGroupPtr->id();
         oldGroupPtr = find(groupId);
         if (oldGroupPtr) {
            returnPtr();
            atomStorage.findGroupAtoms(*oldGroupPtr);
         } else {
            add();
            atomStorage.findGroupAtoms(*newGroupPtr);
         }
      }
      buffer.endRecvBlock();
      assert(buffer.recvSize() == 0);
   }
   #endif // endif ifdef UTIL_MPI

   /*
   * Set ghost communication flags for all atoms in incomplete groups.
   *
   * Precondition: This is called by exchangeAtoms after exchanging atoms 
   * and groups between neighboring processors. At this point, there are
   * no ghosts atoms.
   *
   * Algorithm: Loop over all Group<N> objects in the group storage. 
   * For each group, check if the group is incomplete, implying that one or
   * more atoms in the group are owned by another processor. If the group 
   * is incomplete, loop over 6 transfer directions. For each direction,
   * if the group is marked for sending in that direction, set the ghost
   * ghost communication flag for transfer in that direction for every 
   * local atom in the group. Also add each such atom to sendArray(i, j).
   *
   * Note: If a group is incomplete on this processor, and thus
   * contains atoms owned by other processors, the algorithm assumes
   * that the ghost communication flag for each atom will be set by
   * the processor that owns the atom.
   */
   template <int N> void 
   GroupStorage<N>::markGhosts(AtomStorage& atomStorage, 
                               FMatrix<APArray<Atom>, Dimension, 2>&  sendArray, 
                               IntVector& gridFlags)
   {
      GroupIterator<N> groupIter;
      Atom* atomPtr;
      Plan* planPtr;
      int i, j, k, nAtom;

      // Loop over groups
      begin(groupIter);
      for ( ; groupIter.notEnd(); ++groupIter) {

         #ifdef UTIL_DEBUG
         #ifdef DDMD_GROUP_STORAGE_DEBUG
         // Validate group
         int atomId;
         nAtom  = 0;
         for (k = 0; k < N; ++k) {
            atomPtr = groupIter->atomPtr(k);
            atomId  = groupIter->atomId(k);
            if (atomPtr != 0) {
               if (atomPtr != atomStorage.find(atomId)) {
                  UTIL_THROW("Error in atom pointer in bond");
               }
               if (atomPtr->isGhost()) {
                  UTIL_THROW("Pointer to ghost atom in bond");
               } else {
                  ++nAtom;
               }
            } else { // if atomPtr == 0
               atomPtr = atomStorage.find(atomId);
               if (atomPtr) {
                  if (!atomPtr->isGhost()) {
                     UTIL_THROW("Missing pointer to local atom in bond");
                  }
               }
            }
         }
         assert(nAtom == groupIter->nPtr());
         if (nAtom == 0) {
            UTIL_THROW("Empty group");
         }
         #endif // ifdef DDMD_GROUP_STORAGE_DEBUG
         #endif // ifdef UTIL_DEBUG

         // If this group is incomplete, set ghost flags for atoms 
         nAtom = groupIter->nPtr();
         if (nAtom < N) {
            for (i = 0; i < Dimension; ++i) {
               if (gridFlags[i]) {
                  for (j = 0; j < 2; ++j) {
                     if (groupIter->plan().ghost(i, j)) {
                        for (k = 0; k < N; ++k) {
                           atomPtr = groupIter->atomPtr(k);
                           if (atomPtr) {
                              assert(!atomPtr->isGhost());
                              planPtr = &atomPtr->plan();
                              if (!planPtr->ghost(i, j)) { 
                                 planPtr->setGhost(i, j);
                                 sendArray(i, j).append(*atomPtr);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }

      }
   }

   /*
   * Find ghost members of groups after exchanging all ghosts.
   */
   template <int N>
   void GroupStorage<N>::findGhosts(AtomStorage& atomStorage)
   {
      GroupIterator<N> groupIter;
      int nAtom;
      for (begin(groupIter); groupIter.notEnd(); ++groupIter) {
         nAtom = groupIter->nPtr();
         if (nAtom < N) {
            nAtom = atomStorage.findGroupAtoms(*groupIter);
            if (nAtom < N) {
               UTIL_THROW("Incomplete group after search for ghosts");
            }
         }
      }
   }

} // namespace DdMd
#endif
