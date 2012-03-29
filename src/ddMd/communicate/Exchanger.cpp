#ifndef DDMD_EXCHANGER_CPP
#define DDMD_EXCHANGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Exchanger.h"
#include "Domain.h"
#include "Buffer.h"
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/storage/GroupStorage.h>
#include <ddMd/storage/GroupIterator.h>
#include <util/mpi/MpiLogger.h>
#include <util/global.h>

#include <algorithm>
#include <string>

//#define EXCHANGER_DEBUG

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   Exchanger::Exchanger()
    : pairCutoff_(-1.0)
   {}

   /*
   * Destructor.
   */
   Exchanger::~Exchanger()
   {}

   /*
   * Set pointers to associated objects.
   */
   void Exchanger::associate(const Domain& domain, 
                             const Boundary& boundary, 
                             AtomStorage& atomStorage, 
                             GroupStorage<2>& bondStorage, 
                             #ifdef INTER_ANGLE
                             GroupStorage<3>& angleStorage, 
                             #endif
                             #ifdef INTER_DIHEDRAL
                             GroupStorage<4>& dihedralStorage, 
                             #endif
                             Buffer& buffer)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;
      atomStoragePtr_ = &atomStorage;
      bondStoragePtr_ = &bondStorage;
      #ifdef INTER_ANGLE
      angleStoragePtr_ = &angleStorage;
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_ = &dihedralStorage;
      #endif
      bufferPtr_ = &buffer;
   }

   /*
   * Allocate memory.
   */
   void Exchanger::allocate()
   {
      // Preconditions
      if (!bufferPtr_->isInitialized()) {
         UTIL_THROW("Buffer must be allocated before Exchanger");
      }

      int sendRecvCapacity = bufferPtr_->ghostCapacity()/2;

      // Reserve space for all ghost send and recv arrays
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {
            sendArray_(i, j).reserve(sendRecvCapacity);
            recvArray_(i, j).reserve(sendRecvCapacity);
         }
      }
      emptyBonds_.reserve(sendRecvCapacity);
      #ifdef INTER_ANGLE
      emptyAngles_.reserve(sendRecvCapacity);
      #endif
      #ifdef INTER_DIHEDRAL
      emptyDihedrals_.reserve(sendRecvCapacity);
      #endif
   }

   /*
   * Set slab width used for ghosts.
   */
   void Exchanger::setPairCutoff(double pairCutoff)
   {  pairCutoff_ = pairCutoff; }

   #ifdef UTIL_MPI
   /**
   * Exchange local atoms and ghosts.
   */
   void Exchanger::exchange()
   {
      exchangeAtoms();
      exchangeGhosts();
   }

   template <int N>
   void Exchanger::initGroupGhostPlan (GroupStorage<N>& storage)
   {
      double coordinate;
      GroupIterator<N> groupIter;
      Atom* atomPtr;
      int nIn;
      int nOut;
      int i, j, k;
      bool choose;

      storage.begin(groupIter);
      for ( ; groupIter.notEnd(); ++groupIter) {

         // Compute ghost communication plan for group
         groupIter->plan().clearFlags();
         for (i = 0; i < Dimension; ++i) {
            if (domainPtr_->grid().dimension(i) > 1) {
               for (j = 0; j < 2; ++j) {
                  choose = false;
                  nIn = 0;
                  nOut = 0;
                  for (k = 0; k < N; ++k) {
                     atomPtr = groupIter->atomPtr(k);
                     if (atomPtr) {
                        coordinate = atomPtr->position()[i];
                        if (atomPtr->isGhost()) {
                           if (j == 0) {
                              if (coordinate < inner_(i, j)) {
                                 ++nOut;
                              }
                              if (coordinate > outer_(i, j)) {
                                 ++nIn;
                              }
                           } else {
                              if (coordinate > inner_(i, j)) {
                                 ++nOut;
                              }
                              if (coordinate < outer_(i, j)) {
                                 ++nIn;
                              }
                           }
                        } else { 
                           if (atomPtr->plan().exchange(i, j)) {
                              ++nOut;
                           } else {
                              ++nIn;
                           }
                        }
                     } else {
                        choose = true;
                        break;
                     }
                  } // end for k
                  if (nOut > 0 && nIn > 0) {
                     choose = true;
                  }
                  if (choose) {
                     groupIter->plan().setGhost(i, j);
                  } else {
                     groupIter->plan().clearGhost(i, j);
                  }
               } // end for j
            }
         } // end for i

         // Clear pointers to all ghost atoms in group
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

   /**
   * Pack bonds that contain postmarked atoms.
   */
   template <int N>
   void Exchanger::packGroups(int i, int j, 
                              GroupStorage<N>& storage, 
                              APArray< Group<N> >& emptyGroups)
   {
      GroupIterator<N> groupIter;
      Atom* atomPtr;
      int k, nAtom;
      bool choose;

      bufferPtr_->beginSendBlock(Buffer::GROUP, N);
      emptyGroups.clear();
      storage.begin(groupIter);
      for ( ; !groupIter.atEnd(); ++groupIter) {
         atomStoragePtr_->findGroupAtoms(*groupIter);
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
            emptyGroups.append(*groupIter);
         }
         if (choose) {
            bufferPtr_->packGroup<N>(*groupIter);
         }
      }
      bufferPtr_->endSendBlock();
   }

   /*
   * Remove empty bonds from GroupStorage<N>.
   */
   template <int N>
   void Exchanger::removeEmptyGroups(GroupStorage<N>& storage,
                                     APArray< Group<N> >& emptyGroups)
   {
      Atom* atomPtr;
      int k, m, atomId;
      int nEmpty = emptyGroups.size();
      for (k = 0; k < nEmpty; ++k) {
         #ifdef UTIL_DEBUG
         #ifdef EXCHANGER_DEBUG
         // Confirm that group is actually empty
         for (m = 0; m < N; ++m) {
            atomId  = emptyGroups[k].atomId(m);
            atomPtr = emptyGroups[k].atomPtr(m);
            assert(atomPtr == 0);
            assert(atomStoragePtr_->find(atomId) == 0);
         }
         #endif
         #endif
         storage.remove(&(emptyGroups[k]));
      }
   }

   // Unpack bonds into bondStorage.
   template <int N>
   void Exchanger::unpackGroups(GroupStorage<N>& storage)
   {
      Group<N>* newGroupPtr;
      Group<N>* oldGroupPtr;
      int bondId;

      bufferPtr_->beginRecvBlock();
      while (bufferPtr_->recvSize() > 0) {
         newGroupPtr = storage.newPtr();
         bufferPtr_->unpackGroup<N>(*newGroupPtr);
         bondId = newGroupPtr->id();
         oldGroupPtr = storage.find(bondId);
         if (oldGroupPtr) {
            storage.returnPtr();
            atomStoragePtr_->findGroupAtoms(*oldGroupPtr);
         } else {
            storage.add();
            atomStoragePtr_->findGroupAtoms(*newGroupPtr);
         }
      }
      assert(bufferPtr_->recvSize() == 0);
   }

   template <int N>
   void Exchanger::finishGroupGhostPlan(GroupStorage<N>& storage)
   {
      GroupIterator<N> groupIter;
      Atom* atomPtr;
      int i, j, k, atomId, nAtom;
      bool choose;

      // Set ghost communication flags for atoms in incomplete bonds
      storage.begin(groupIter);
      for ( ; groupIter.notEnd(); ++groupIter) {

         #ifdef UTIL_DEBUG
         #ifdef EXCHANGER_DEBUG
         nAtom  = 0;
         for (k = 0; k < N; ++k) {
            atomPtr = groupIter->atomPtr(k);
            atomId  = groupIter->atomId(k);
            if (atomPtr != 0) {
               if (atomPtr != atomStoragePtr_->find(atomId)) {
                  UTIL_THROW("Error in atom pointer in bond");
               }
               if (atomPtr->isGhost()) {
                  UTIL_THROW("Pointer to ghost atom in bond");
               } else {
                  ++nAtom;
               }
            } else { // if atomPtr == 0
               atomPtr = atomStoragePtr_->find(atomId);
               if (atomPtr) {
                  if (!atomPtr->isGhost()) {
                     UTIL_THROW("Missing pointer to local atom in bond");
                  }
               }
            }
         }
         assert(nAtom == groupIter->nPtr());
         if (nAtom == 0) {
            UTIL_THROW("Empty bond");
         }
         #endif
         #endif

         // Set communication flags for atoms in incomplete groups
         nAtom = groupIter->nPtr();
         if (nAtom < N) {
            for (i = 0; i < Dimension; ++i) {
               if (domainPtr_->grid().dimension(i) > 1) {
                  for (j = 0; j < 2; ++j) {
                     // choose = true;
                     choose = groupIter->plan().ghost(i, j);
                     if (choose) {
                        for (k = 0; k < N; ++k) {
                           atomPtr = groupIter->atomPtr(k);
                           if (atomPtr) {
                              if (!atomPtr->isGhost()) {
                                 atomPtr->plan().setGhost(i, j);
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

   /**
   * Exchange ownership of local atoms.
   *
   * This method should be called before rebuilding the neighbor list on
   * each processor, to exchange ownership of local atoms.
   */
   void Exchanger::exchangeAtoms()
   {
      Vector lengths = boundaryPtr_->lengths();
      double bound, inner, outer, middle;
      double coordinate;
      AtomIterator atomIter;
      GhostIterator ghostIter;
      Atom* atomPtr;
      int i, j, jc, k, m, source, dest, nSend;
      int atomId, nAtom;
      int myRank = domainPtr_->gridRank();
      int shift;
      bool choose;

      // Set domain and slab boundaries
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {
            // j = 0 sends to lower coordinate i, bound is minimum
            // j = 1 sends to higher coordinate i, bound is maximum
            bound = domainPtr_->domainBound(i, j);
            bound_(i, j) = bound;
            if (j == 0) { // Communicate with lower index
               inner_(i,j) = bound + pairCutoff_;
               outer_(i, j)= bound - pairCutoff_;
            } else { // j == 1, communicate with upper index
               inner_(i, j) = bound - pairCutoff_;
               outer_(i, j) = bound + pairCutoff_;
            }
         }
      }

      // Compute communication plans for all local atoms
      atomStoragePtr_->begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {

         atomIter->plan().clearFlags();

         // Cartesian directions
         for (i = 0; i < Dimension; ++i) {
  
            coordinate = atomIter->position()[i];
 
            // Transmission direction
            for (j = 0; j < 2; ++j) {
   
               // j = 0 sends to lower coordinate i
               // j = 1 sends to higher coordinate i

               // Index for conjugate (reverse) direction
               if (j == 0) jc = 1;
               if (j == 1) jc = 0;
   
               if (j == 0) { // Communicate with lower index
                  if (coordinate < bound_(i, j)) {
                     atomIter->plan().setExchange(i, j);
                     if (coordinate > outer_(i, j)) {
                        atomIter->plan().setGhost(i, jc);
                     }
                  } else {
                     if (coordinate < inner_(i, j)) {
                        atomIter->plan().setGhost(i, j);
                     }
                  }
               } else { // j == 1, communicate with upper index
                  if (coordinate > bound_(i, j)) {
                     atomIter->plan().setExchange(i, j);
                     if (coordinate < outer_(i, j)) {
                        atomIter->plan().setGhost(i, jc);
                     }
                  } else {
                     if (coordinate > inner_(i, j)) {
                        atomIter->plan().setGhost(i, j);
                     }
                  }
               }
   
            } // end for j
         } // end for i

      } // end atom loop, end compute plan

      initGroupGhostPlan<2>(*bondStoragePtr_);
      #ifdef INTER_ANGLE
      initGroupGhostPlan<3>(*angleStoragePtr_);
      #endif
      #ifdef INTER_DIHEDRAL
      initGroupGhostPlan<4>(*dihedralStoragePtr_);
      #endif

      // Clear all ghost atoms from AtomStorage
      atomStoragePtr_->clearGhosts();

      #ifdef UTIL_DEBUG
      #ifdef EXCHANGER_DEBUG
      int nAtomTotal;
      atomStoragePtr_->computeNAtomTotal(domainPtr_->communicator());
      if (myRank == 0) {
         nAtomTotal = atomStoragePtr_->nAtomTotal();
      }
      bondStoragePtr_->isValid(*atomStoragePtr_, 
                               domainPtr_->communicator(), false);
      #endif
      #endif

      // Cartesian directions
      for (i = 0; i < Dimension; ++i) {

         // Transmission direction
         for (j = 0; j < 2; ++j) {

            // j = 0 sends to lower  coordinate i
            // j = 1 sends to higher coordinate i

            // Index for conjugate (reverse) direction
            if (j == 0) jc = 1;
            if (j == 1) jc = 0;

            source = domainPtr_->sourceRank(i, j); // rank to receive from
            dest = domainPtr_->destRank(i, j);     // rank to send to
            bound = domainPtr_->domainBound(i, j); // bound for send
            inner = inner_(i, jc);             // inner bound upon receipt 
            shift = domainPtr_->shift(i, j);   // shift for periodic b.c.

            sendArray_(i, j).clear();

            if (domainPtr_->grid().dimension(i) > 1) {
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::ATOM);
            }

            // Choose atoms for sending, pack and mark for removal.
            atomStoragePtr_->begin(atomIter);
            for ( ; !atomIter.atEnd(); ++atomIter) {

               #ifdef UTIL_DEBUG
               #ifdef EXCHANGER_DEBUG
               if (j == 0) {
                  choose = (atomIter->position()[i] < bound);
               } else {
                  choose = (atomIter->position()[i] > bound);
               }
               assert(choose == atomIter->plan().exchange(i, j));
               #endif
               #endif

               choose = atomIter->plan().exchange(i, j);
               if (choose) {

                  if (domainPtr_->grid().dimension(i) > 1) {

                     sendArray_(i, j).append(*atomIter);
                     bufferPtr_->packAtom(*atomIter);

                  } else {

                     // Shift position if required by periodic b.c.
                     if (shift) {
                        atomIter->position()[i] += shift * lengths[i];
                     }
                     assert(atomIter->position()[i] 
                            > domainPtr_->domainBound(i, 0));
                     assert(atomIter->position()[i] 
                            < domainPtr_->domainBound(i, 1));

                     // For gridDimension==1, only nonbonded ghosts exist.
                     // The following assertion applies to these.
                     assert(!atomIter->plan().ghost(i, j));

                     #if UTIL_DEBUG
                     // Check ghost communication plan
                     if (j == 0 && atomIter->position()[i] > inner) { 
                        assert(atomIter->plan().ghost(i, 1));
                     } else 
                     if (j == 1 && atomIter->position()[i] < inner) {
                        assert(atomIter->plan().ghost(i, 0));
                     }
                     #endif

                  }
               }

            } // end atom loop

            // Send and receive only if dimension(i) > 1
            if (domainPtr_->grid().dimension(i) > 1) {

               // End atom send block
               bufferPtr_->endSendBlock();

               // Pack bonds that contain postmarked atoms.
               packGroups<2>(i, j, *bondStoragePtr_, emptyBonds_);
               #ifdef INTER_ANGLE
               packGroups<3>(i, j, *angleStoragePtr_, emptyAngles_);
               #endif
               #ifdef INTER_DIHEDRAL
               packGroups<4>(i, j, *dihedralStoragePtr_, emptyDihedrals_);
               #endif

               /*
               * Note: Removal cannot be done within the above loops over 
               * atoms and groups because element removal invalidates a
               * PArray iterator.
               */

               // Remove chosen atoms from atomStorage
               nSend = sendArray_(i, j).size();
               for (k = 0; k < nSend; ++k) {
                  atomStoragePtr_->removeAtom(&sendArray_(i, j)[k]);
               }
     
               // Remove empty bonds from bondStorage.
               removeEmptyGroups<2>(*bondStoragePtr_, emptyBonds_);
               #ifdef INTER_ANGLE
               removeEmptyGroups<3>(*angleStoragePtr_, emptyAngles_);
               #endif
               #ifdef INTER_DIHEDRAL
               removeEmptyGroups<4>(*dihedralStoragePtr_, emptyDihedrals_);
               #endif

               // Send to processor dest and receive from processor source
               bufferPtr_->sendRecv(domainPtr_->communicator(), 
                                    source, dest);

               // Unpack atoms into atomStorage
               bufferPtr_->beginRecvBlock();
               while (bufferPtr_->recvSize() > 0) {
                  atomPtr = atomStoragePtr_->newAtomPtr();
                  bufferPtr_->unpackAtom(*atomPtr);

                  if (shift) {
                     atomPtr->position()[i] += shift * lengths[i];
                  }
                  assert(atomPtr->position()[i] 
                         > domainPtr_->domainBound(i, 0));
                  assert(atomPtr->position()[i] 
                         < domainPtr_->domainBound(i, 1));


                  #ifdef UTIL_DEBUG
                  // Check ghost plan
                  assert(!atomPtr->plan().ghost(i, j));
                  if (j == 0) {
                     assert( atomPtr->plan().ghost(i, 1)
                             == (atomPtr->position()[i] > inner) );
                  } else {
                     assert( atomPtr->plan().ghost(i, 0)
                             == (atomPtr->position()[i] < inner) );
                  }
                  #endif

                  atomStoragePtr_->addNewAtom();
               }
               assert(bufferPtr_->recvSize() == 0);

               // Unpack bonds into bondStorage.
               unpackGroups<2>(*bondStoragePtr_);
               #ifdef INTER_ANGLE
               unpackGroups<3>(*angleStoragePtr_);
               #endif
               #ifdef INTER_DIHEDRAL
               unpackGroups<4>(*dihedralStoragePtr_);
               #endif

            } // end if gridDimension > 1

         } // end for j (direction 0, 1)
      } // end for i (Cartesian index)

      /*
      * At this point:
      * No ghost atoms exist.
      * All atoms are on correct processor.
      * No Groups are empty.
      * All pointer to local atoms in Groups are set.
      * All pointers to ghost atoms in Groups are null.
      */

      #ifdef UTIL_DEBUG
      #ifdef EXCHANGER_DEBUG
      atomStoragePtr_->computeNAtomTotal(domainPtr_->communicator());
      if (myRank == 0) {
         assert(nAtomTotal = atomStoragePtr_->nAtomTotal());
      }
      atomStoragePtr_->isValid();
      bondStoragePtr_->isValid(*atomStoragePtr_, 
                               domainPtr_->communicator(), false);
      #ifdef INTER_ANGLE
      angleStoragePtr_->isValid(*atomStoragePtr_, 
                                domainPtr_->communicator(), false);
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_->isValid(*atomStoragePtr_, 
                                   domainPtr_->communicator(), false);
      #endif
      #endif
      #endif

      // Set ghost communication flags for atoms in incomplete bonds
      finishGroupGhostPlan<2>(*bondStoragePtr_);
      #ifdef INTER_ANGLE
      finishGroupGhostPlan<3>(*angleStoragePtr_);
      #endif
      #ifdef INTER_DIHEDRAL
      finishGroupGhostPlan<4>(*dihedralStoragePtr_);
      #endif

   }

   /*
   * Find all ghost members of groups at the end of exchangeGhosts.
   */
   template <int N>
   void Exchanger::findGroupGhosts(GroupStorage<N>& storage)
   {
      GroupIterator<N> groupIter;
      int nAtom;
      storage.begin(groupIter);
      for ( ; !groupIter.atEnd(); ++groupIter) {
         nAtom = groupIter->nPtr();
         if (nAtom < N) {
            nAtom = atomStoragePtr_->findGroupAtoms(*groupIter);
            if (nAtom < N) {
               UTIL_THROW("Incomplete group after search for ghosts");
            }
         }
      }
   }

   /*
   * Exchange ghost atoms.
   *
   * Call immediately after exchangeAtoms and before rebuilding the 
   * neighbor list on time steps that require reneighboring. Uses
   * ghost communication plans computed in exchangeAtoms.
   */
   void Exchanger::exchangeGhosts()
   {

      // Preconditions
      assert(bufferPtr_->isInitialized());
      assert(domainPtr_->isInitialized());
      assert(domainPtr_->hasBoundary());

      Vector        lengths = boundaryPtr_->lengths();
      double        bound, inner, coord;
      AtomIterator  localIter;
      GhostIterator ghostIter;
      Atom* atomPtr;
      int i, j, jc, source, dest, shift;
      int myRank = domainPtr_->gridRank();

      // Check that all ghosts are cleared upon entry.
      if (atomStoragePtr_->nGhost() != 0) {
         UTIL_THROW("atomStoragePtr_->nGhost() != 0");
      }

      // Cartesian directions
      for (i = 0; i < Dimension; ++i) {

         // Transmit directions
         for (j = 0; j < 2; ++j) {

            // j = 0: Send ghosts near minimum boundary to lower coordinate
            // j = 1: Send ghosts near maximum boundary to higher coordinate

            // Set index for reverse direction
            if (j == 0) jc = 1;
            if (j == 1) jc = 0;

            bound = bound_(i, j);
            inner = inner_(i, j);
       
            // Shift on receiving processor for periodic boundary conditions
            shift = domainPtr_->shift(i, j);

            if (domainPtr_->grid().dimension(i) > 1)  {
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::GHOST);
            }

            sendArray_(i, j).clear();
            recvArray_(i, j).clear();

            // Loop over local atoms on this processor
            atomStoragePtr_->begin(localIter);
            for ( ; !localIter.atEnd(); ++localIter) {

               coord = localIter->position()[i];

               #if 0 
               // This test is valid only for nonbonded ghosts.
               {
                  bool choose;
                  if (j == 0) {
                     assert(coord > bound);
                     choose = (coord < inner);
                  } else {
                     assert(coord < bound);
                     choose = (coord > inner);
                  }
                  assert(choose == localIter->plan().ghost(i, j));
               }
               #endif

               // If this local atom is selected for sending as a ghost.
               if ( localIter->plan().ghost(i, j) ) {

                  sendArray_(i, j).append(*localIter);

                  if (domainPtr_->grid().dimension(i) > 1)  {

                     // Pack atom for sending 
                     bufferPtr_->packGhost(*localIter);

                  } else {  // if grid dimension == 1

                     // Make a ghost copy of local atom on this processor
                     atomPtr = atomStoragePtr_->newGhostPtr();
                     recvArray_(i, j).append(*atomPtr);
                     atomPtr->setId(localIter->id());
                     atomPtr->setTypeId(localIter->typeId());
                     atomPtr->plan().setFlags(localIter->plan().flags());
                     atomPtr->position() = localIter->position();
                     if (shift) {
                        atomPtr->position()[i] += shift * lengths[i];
                     }
                     atomStoragePtr_->addNewGhost();

                     // Prevent ghost copy from being re-copied within
                     // following loop over ghosts on same processor.
                     atomPtr->plan().clearGhost(i, j);

                     #ifdef UTIL_DEBUG
                     // Validate shifted ghost coordinate 
                     if (j == 0) {
                        assert(atomPtr->position()[i] > bound_(i, 1));
                     } else {
                        assert(atomPtr->position()[i] < bound_(i, 0));
                     }
                     #endif

                  }

               }

            }

            // Loop over ghosts on this processor, for resending.
            atomStoragePtr_->begin(ghostIter);
            for ( ; !ghostIter.atEnd(); ++ghostIter) {

               coord = ghostIter->position()[i];

               #if 0
               // This test is valid only for nonbonded ghosts.
               {
                  bool choose;
                  if (j == 0) {
                     choose = (coord > bound) && (coord < inner);
                  } else {
                     choose = (coord < bound) && (coord > inner);
                  }
                  assert(choose == ghostIter->plan().ghost(i, j));
               }
               #endif

               // If this ghost is selected for re-sending
               if ( ghostIter->plan().ghost(i, j) ) {

                  sendArray_(i, j).append(*ghostIter);

                  if (domainPtr_->grid().dimension(i) > 1)  {
                     
                     // Pack ghost for resending
                     bufferPtr_->packGhost(*ghostIter);

                  } else {  // if grid dimension == 1

                     // Make another ghost copy on the same processor
                     atomPtr = atomStoragePtr_->newGhostPtr();
                     recvArray_(i, j).append(*atomPtr);
                     atomPtr->setId(ghostIter->id());
                     atomPtr->setTypeId(ghostIter->typeId());
                     atomPtr->plan().setFlags(ghostIter->plan().flags());
                     atomPtr->position() = ghostIter->position();
                     if (shift) {
                        atomPtr->position()[i] += shift * lengths[i];
                     }
                     atomStoragePtr_->addNewGhost();

                     #ifdef EXCHANGE_DEBUG
                     #ifdef UTIL_DEBUG
                     // Validate shifted position
                     if (j == 0) {
                        assert(atomPtr->position()[i] > bound_(i, 1));
                     } else {
                        assert(atomPtr->position()[i] < bound_(i, 0));
                     }
                     #endif
                     #endif

                  }

               }

            }

            // Send and receive buffers
            if (domainPtr_->grid().dimension(i) > 1) {

               bufferPtr_->endSendBlock();

               source = domainPtr_->sourceRank(i, j);
               dest   = domainPtr_->destRank(i, j);
               bufferPtr_->sendRecv(domainPtr_->communicator(), source, dest);

               // Unpack ghosts and add to recvArray
               bufferPtr_->beginRecvBlock();
               while (bufferPtr_->recvSize() > 0) {

                  atomPtr = atomStoragePtr_->newGhostPtr();
                  bufferPtr_->unpackGhost(*atomPtr);
                  if (shift) {
                     atomPtr->position()[i] += shift * lengths[i];
                  }
                  recvArray_(i, j).append(*atomPtr);
                  atomStoragePtr_->addNewGhost();

                  // Prohibit sending back ghost in reverse direction
                  if (j == 0) {
                     atomPtr->plan().clearGhost(i, 1);
                  }

                  #ifdef UTIL_DEBUG
                  // Validate ghost coordinate on receiving processor.
                  if (j == 0) {
                     if (atomPtr->position()[i] < bound_(i, 1)) {
                        std::cout << i << "  " << j 
                                  << "  " << atomPtr->plan()
                                  << "  " << atomPtr->position()[i]
                                  << "  " << bound_(i, 0)
                                  << std::endl;
                     }
                     //assert(atomPtr->position()[i] > bound_(i, 1));
                  } else {
                     if (atomPtr->position()[i] > bound_(i, 0)) {
                        std::cout << i << "  " << j 
                                  << "  " << atomPtr->plan()
                                  << "  " << atomPtr->position()[i]
                                  << "  " << bound_(i, 0)
                                  << std::endl;
                     }
                     //assert(atomPtr->position()[i] < bound_(i, 0));
                  }
                  #endif

               }

            }

         } // end for transmit direction j = 0, 1

      } // end for Cartesian index i


      // Find ghost atoms for all incomplete bonds
      findGroupGhosts<2>(*bondStoragePtr_);
      #ifdef INTER_ANGLE
      findGroupGhosts<3>(*angleStoragePtr_);
      #endif
      #ifdef INTER_DIHEDRAL
      findGroupGhosts<4>(*dihedralStoragePtr_);
      #endif

      #ifdef UTIL_DEBUG
      #ifdef EXCHANGER_DEBUG
      atomStoragePtr_->isValid();
      bondStoragePtr_->isValid(*atomStoragePtr_, domainPtr_->communicator(),
                               true);
      #endif
      #endif

   }

   /*
   * Update ghost atom coordinates.
   *
   * Call on time steps for which no reneighboring is required. 
   */
   void Exchanger::update()
   {
      Vector lengths = boundaryPtr_->lengths();
      Atom*  atomPtr;
      int    i, j, k, source, dest, size, shift;

      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {

            // Shift on receiving processor for periodic boundary conditions
            shift = domainPtr_->shift(i, j);

            if (domainPtr_->grid().dimension(i) > 1) {

               // Pack ghost positions for sending
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::UPDATE);
               size = sendArray_(i, j).size();
               for (k = 0; k < size; ++k) {
                  bufferPtr_->packUpdate(sendArray_(i, j)[k]);
               }
               bufferPtr_->endSendBlock();
  
               // Send and receive buffers
               source = domainPtr_->sourceRank(i, j);
               dest   = domainPtr_->destRank(i, j);
               bufferPtr_->sendRecv(domainPtr_->communicator(), source, dest);
   
               // Unpack ghost positions
               bufferPtr_->beginRecvBlock();
               size = recvArray_(i, j).size();
               for (k = 0; k < size; ++k) {
                  atomPtr = &recvArray_(i, j)[k];
                  bufferPtr_->unpackUpdate(*atomPtr);
                  if (shift) {
                     atomPtr->position()[i] += shift * lengths[i];
                  }
               }

            } else {

               // If grid().dimension(i) == 1, then copy positions of atoms
               // listed in sendArray to those listed in the recvArray.

               size = sendArray_(i, j).size();
               assert(size == recvArray_(i, j).size());
               for (k = 0; k < size; ++k) {
                  atomPtr = &recvArray_(i, j)[k];
                  atomPtr->position() = sendArray_(i, j)[k].position();
                  if (shift) {
                     atomPtr->position()[i] += shift * lengths[i];
                  }
               }

            }  

         } // transmit direction j = 0, 1

      } // Cartesian direction i = 0, ..., Dimension - 1

   }
   #endif

}
#endif
