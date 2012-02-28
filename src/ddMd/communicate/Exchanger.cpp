#ifndef EXCHANGER_CPP
#define EXCHANGER_CPP

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
#include <ddMd/storage/BondStorage.h>
#include <ddMd/storage/GroupIterator.h>
#include <util/mpi/MpiLogger.h>
#include <util/global.h>

#include <algorithm>
#include <string>

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
                             BondStorage& bondStorage, 
                             Buffer& buffer)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;
      atomStoragePtr_ = &atomStorage;
      bondStoragePtr_ = &bondStorage;
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
      incompleteBonds_.reserve(sendRecvCapacity);
   }

   /*
   * Set slab width used for ghosts.
   */
   void Exchanger::setPairCutoff(double pairCutoff)
   {  pairCutoff_ = pairCutoff; }

   #ifdef UTIL_MPI
   /**
   * Exchange ownership of local atoms.
   *
   * This method should be called before rebuilding the neighbor list on
   * each processor, to exchange ownership of local atoms.
   */
   void Exchanger::exchangeAtoms()
   {
      Vector lengths = boundaryPtr_->lengths();
      double bound, inner, coordinate;
      AtomIterator atomIter;
      GhostIterator ghostIter;
      GroupIterator<2> bondIter;
      Atom* atomPtr;
      int i, j, jc, k, m, source, dest, nSend;
      int atomId, nAtom;
      int myRank = domainPtr_->gridRank();
      int shift;
      bool choose;

      #if 1
      // Calculate atom communication plans for all local atoms
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
   
               // Boundary (j=0 -> minimum, j=1 -> maximum)
               bound = domainPtr_->domainBound(i, j);

               // Decide upon plan  direction i, j
               if (j == 0) { // Communicate with lower index
                  inner = bound + pairCutoff_;
                  if (coordinate < inner) {
                     atomIter->plan().setGhost(i, j);
                     if (coordinate < bound) {
                        atomIter->plan().setExchange(i, j);
                     }
                  }
               } else { // j == 1, communicate with upper index
                  inner = bound - pairCutoff_;
                  if (coordinate > inner) {
                     atomIter->plan().setGhost(i, j);
                     if (coordinate > bound) {
                        atomIter->plan().setExchange(i, j);
                     }
                  }
               }
   
            } // end loop j
         } // loop i

      } // end atom loop, end compute plan
      #endif

      // In all bonds, clear pointers to ghosts and check local pointers.
      bondStoragePtr_->begin(bondIter);
      for ( ; bondIter.notEnd(); ++bondIter) {
         nAtom = 0;
         for (k = 0; k < 2; ++k) {
            atomPtr = bondIter->atomPtr(k);
            atomId = bondIter->atomId(k);
            if (atomPtr != 0) {
               #ifdef UTIL_DEBUG
               if (atomPtr != atomStoragePtr_->find(atomId)) {
                  UTIL_THROW("Error in atom pointer in bond");
               }
               #endif
               if (atomPtr->isGhost()) {
                  bondIter->clearAtomPtr(k);
               } else {
                  atomId = atomPtr->id();
                  ++nAtom;
               }
            }
            #ifdef UTIL_DEBUG 
            else { // if atomPtr == 0
               atomPtr = atomStoragePtr_->find(atomId);
               if (atomPtr) {
                  if (!atomPtr->isGhost()) {
                     UTIL_THROW("Missing pointer to local atom in bond");
                  }
               }
            }
            #endif
         }
         if (nAtom == 0) {
            UTIL_THROW("Empty bond");
         }
      }

      // Clear ghost atoms from AtomStorage
      atomStoragePtr_->clearGhosts();

      #if 0
      int nAtomTotal; // Total number of atoms on all processors
      if (myRank == 0) {
         nAtomTotal = atomStoragePtr_->nAtomTotal();
      }
      bondStoragePtr_->isValid(*atomStoragePtr_, domainPtr_->communicator(), 
                               false);
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

            //Processor to receive from
            source = domainPtr_->sourceRank(i, j);

            //Processor to send to
            dest = domainPtr_->destRank(i, j);

            // Boundary for sending processor (j=0 -> minimum, j=1 -> maximum)
            bound = domainPtr_->domainBound(i, j);

            // Inner slab boundary for receiving processor. 
            if (j == 0) {
               inner = domainPtr_->domainBound(i, 1) - pairCutoff_;
            } else {
               inner = domainPtr_->domainBound(i, 0) + pairCutoff_;
            }

            // Shift due to periodic boundary conditions
            shift = domainPtr_->shift(i, j);

            sendArray_(i, j).clear();

            if (domainPtr_->grid().dimension(i) > 1) {
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::ATOM);
            }

            // Choose atoms for sending, pack and mark for removal.
            atomStoragePtr_->begin(atomIter);
            for ( ; !atomIter.atEnd(); ++atomIter) {

               //atomIter->plan().clearExchange(i, j);

               if (j == 0) {
                  choose = (atomIter->position()[i] < bound);
               } else {
                  choose = (atomIter->position()[i] > bound);
               }
               assert(choose == atomIter->plan().exchange(i, j));

               if (choose) {

                  if (domainPtr_->grid().dimension(i) > 1) {

                     sendArray_(i, j).append(*atomIter);
                     bufferPtr_->packAtom(*atomIter);
                     atomIter->plan().setExchange(i, j);

                  } else {

                     // Shift position if required by periodic b.c.
                     if (shift) {
                        atomIter->position()[i] += shift * lengths[i];
                     }
                     if (atomIter->plan().ghost(i, j)) {
                        atomIter->plan().clearGhost(i, j);
                     }
                     if (j == 0 && atomIter->position()[i] > inner) { 
                           atomIter->plan().setGhost(i, 1);
                     } else 
                     if (j == 1 && atomIter->position()[i] < inner) { 
                           atomIter->plan().setGhost(i, 0);
                     }
                     assert(atomIter->position()[i] 
                            > domainPtr_->domainBound(i, 0));
                     assert(atomIter->position()[i] 
                            < domainPtr_->domainBound(i, 1));

                  }
               }

            } // end atom loop

            // Send and receive only if dimension(i) > 1
            if (domainPtr_->grid().dimension(i) > 1) {

               // End atom send block
               bufferPtr_->endSendBlock();

               // Loop over bonds
               // Pack bonds with postmarked atoms.
               bufferPtr_->beginSendBlock(Buffer::GROUP, 2);
               emptyBonds_.clear();
               bondStoragePtr_->begin(bondIter);
               for ( ; !bondIter.atEnd(); ++bondIter) {
                  atomStoragePtr_->findGroupAtoms(*bondIter);
                  choose = false;
                  nAtom = 0;
                  for (k = 0; k < 2; ++k) {
                     atomPtr = bondIter->atomPtr(k);
                     if (atomPtr) {
                        if (atomPtr->plan().exchange(i, j)) {
                           choose = true;
                           bondIter->clearAtomPtr(k);
                        } else {
                           ++nAtom;
                        }
                     }
                  }
                  if (nAtom == 0) {
                     emptyBonds_.append(*bondIter);
                  }
                  if (choose) {
                     bufferPtr_->packGroup<2>(*bondIter);
                  }
               }
               bufferPtr_->endSendBlock();

               /*
               * Note: Removal cannot be done within the above loops over 
               * atoms and groups because element removal invalidates any
               * PArray iterator.
               */

               // Remove chosen atoms from atomStorage
               nSend = sendArray_(i, j).size();
               for (k = 0; k < nSend; ++k) {
                  atomStoragePtr_->removeAtom(&sendArray_(i, j)[k]);
               }

               // Remove empty bonds from bondStorage.
               int nEmpty = emptyBonds_.size();
               for (k = 0; k < nEmpty; ++k) {
                  #if 0
                  #ifdef UTIL_DEBUG
                  // Confirm that bond is empty
                  for (m = 0; m < 2; ++m) {
                     atomId  = emptyBonds_[k].atomId(m);
                     atomPtr = emptyBonds_[k].atomPtr(m);
                     assert(atomPtr == 0);
                     assert(atomStoragePtr_->find(atomId) == 0);
                  }
                  #endif
                  #endif
                  bondStoragePtr_->remove(&(emptyBonds_[k]));
               }

               // Send to processor dest and receive from processor source
               bufferPtr_->sendRecv(domainPtr_->communicator(), source, dest);

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

                  // Adjust ghost plan
                  if (atomPtr->plan().ghost(i, j)) {
                     atomPtr->plan().clearGhost(i, j);
                  }
                  if (j == 0 && atomPtr->position()[i] > inner) { 
                        atomPtr->plan().setGhost(i, 1);
                  } else 
                  if (j == 1 && atomPtr->position()[i] < inner) { 
                        atomPtr->plan().setGhost(i, 0);
                  }

                  atomStoragePtr_->addNewAtom();
               }
               assert(bufferPtr_->recvSize() == 0);

               // Unpack bonds into bondStorage.
               int bondId;
               Group<2>* newBondPtr;
               Group<2>* oldBondPtr;
               bufferPtr_->beginRecvBlock();
               while (bufferPtr_->recvSize() > 0) {
                  newBondPtr = bondStoragePtr_->newPtr();
                  bufferPtr_->unpackGroup<2>(*newBondPtr);
                  bondId = newBondPtr->id();
                  oldBondPtr = bondStoragePtr_->find(bondId);
                  if (oldBondPtr) {
                     bondStoragePtr_->returnPtr();
                     atomStoragePtr_->findGroupAtoms(*oldBondPtr);
                  } else {
                     bondStoragePtr_->add();
                     atomStoragePtr_->findGroupAtoms(*newBondPtr);
                  }
               }
               assert(bufferPtr_->recvSize() == 0);

            } // end if gridDimension > 1

         } // end for j (direction 0, 1)
      } // end for i (Cartesian index)

      #if 0
      atomStoragePtr_->computeNAtomTotal(domainPtr_->communicator());
      if (myRank == 0) {
         assert(nAtomTotal = atomStoragePtr_->nAtomTotal());
      }
      #endif

      #if 0
      // Find atoms for all bonds.
      bondStoragePtr_->begin(bondIter); 
      for ( ; bondIter.notEnd(); ++bondIter) {
         nAtom = atomStoragePtr_->findGroupAtoms(*bondIter);
         assert(nAtom > 0);
      }
      #endif

      #if 0
      bondStoragePtr_->isValid(*atomStoragePtr_, domainPtr_->communicator(),
                               false);
      #endif

      /*
      * At this point:
      * No ghost atoms exist.
      * All atoms are on correct processor.
      * All pointer to local atoms in Groups are set.
      * All pointers to ghost atoms in Groups are null.
      */

   }

   /*
   * Exchange ghost atoms.
   *
   * Call after exchangeAtoms and before rebuilding the neighbor list on
   * time steps that require reneighboring.
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
      int i, j, source, dest, shift;
      int myRank = domainPtr_->gridRank();
      bool choose;

      // Check that all ghosts are cleared upon entry.
      if (atomStoragePtr_->nGhost() > 0) {
         UTIL_THROW("atomStoragePtr_->nGhost() > 0");
      }

      // Cartesian directions
      for (i = 0; i < Dimension; ++i) {

         // Transmit directions
         for (j = 0; j < 2; ++j) {

            // j = 0: Send ghosts near minimum boundary to lower coordinate
            // j = 1: Send ghosts near maximum boundary to higher coordinate

            if (j == 0) {
               bound = domainPtr_->domainBound(i, 0);
               inner = bound + pairCutoff_;
            } else {
               bound = domainPtr_->domainBound(i, 1);
               inner = bound - pairCutoff_;
            }

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

               // Decide whether to send as ghost.
               coord = localIter->position()[i];
               if (j == 0) {
                  assert(coord > bound);
                  choose = (coord < inner);
               } else {
                  assert(coord < bound);
                  choose = (coord > inner);
               }

               //assert(choose == localIter->plan().ghost(i, j));

               #if 0
               if (choose != localIter->plan().ghost(i, j)) {
                  std::cout << "Proc " << myRank << "  "
                            << "atom " << localIter->id() << "  "
                            << "Dir  " << i << "  " << j << "  "
                            << localIter->position() << "  ";
                  for (int a= 0; a < Dimension; ++a) {
                     for (int b = 0; b < 2;  ++b) {
                        std::cout << localIter->plan().exchange(a, b);
                     }
                  }
                  std::cout << "  ";
                  for (int a= 0; a < Dimension; ++a) {
                     for (int b = 0; b < 2;  ++b) {
                        std::cout << localIter->plan().ghost(a, b);
                     }
                  }
                  std::cout << "  " << choose << "  " 
                            << localIter->plan().ghost(i, j);
                  std::cout << std::endl;
                  
                  UTIL_THROW("Assert failed in exchange ghosts");
               }
               #endif

               if (choose) {

                  sendArray_(i, j).append(*localIter);

                  if (domainPtr_->grid().dimension(i) > 1)  {

                     // Pack atom for sending 
                     bufferPtr_->packGhost(*localIter);

                  } else {  // if grid dimension == 1

                     // Make a ghost copy of the local atom on this processor
                     atomPtr = atomStoragePtr_->newGhostPtr();
                     recvArray_(i, j).append(*atomPtr);
                     atomPtr->position() = localIter->position();
                     atomPtr->setTypeId(localIter->typeId());
                     atomPtr->setId(localIter->id());
                     if (shift) {
                        atomPtr->position()[i] += shift * lengths[i];
                     }
                     atomStoragePtr_->addNewGhost();

                     #ifdef UTIL_DEBUG
                     // Validate shifted positions
                     if (j == 0) {
                        assert(atomPtr->position()[i] 
                               > domainPtr_->domainBound(i, 1));
                     } else {
                        assert(atomPtr->position()[i] 
                               < domainPtr_->domainBound(i, 0));
                     }
                     #endif

                  }

               } 

            }

            // Loop over ghosts on this processor, for resending.
            atomStoragePtr_->begin(ghostIter);
            for ( ; !ghostIter.atEnd(); ++ghostIter) {

               // Decide whether to resend this ghost
               coord = ghostIter->position()[i];
               if (j == 0) {
                  choose = (coord > bound) && (coord < inner);
               } else {
                  choose = (coord < bound) && (coord > inner);
               }

               if (choose) {

                  sendArray_(i, j).append(*ghostIter);

                  if (domainPtr_->grid().dimension(i) > 1)  {
                     
                     // Pack ghost for resending
                     bufferPtr_->packGhost(*ghostIter);

                  } else {  // if grid dimension == 1

                     // Make another ghost copy on the same processor
                     atomPtr = atomStoragePtr_->newGhostPtr();
                     recvArray_(i, j).append(*atomPtr);
                     atomPtr->position() = ghostIter->position();
                     atomPtr->setTypeId(ghostIter->typeId());
                     atomPtr->setId(ghostIter->id());
                     if (shift) {
                        atomPtr->position()[i] += shift * lengths[i];
                     }
                     atomStoragePtr_->addNewGhost();

                     #ifdef UTIL_DEBUG
                     // Validate shifted position
                     if (j == 0) {
                        assert(atomPtr->position()[i] 
                               > domainPtr_->domainBound(i, 1));
                     } else {
                        assert(atomPtr->position()[i] 
                               < domainPtr_->domainBound(i, 0));
                     }
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

                  #ifdef UTIL_DEBUG
                  // Validate ghost coordinates on receiving processor.
                  if (j == 0) {
                     assert(atomPtr->position()[i] 
                            > domainPtr_->domainBound(i, 1));
                  } else {
                     assert(atomPtr->position()[i] 
                            < domainPtr_->domainBound(i, 0));
                  }
                  #endif

               }

            }

         } // transmit direction j = 0, 1

      } // Cartesian index i = 0, ..., Dimension - 1

   }

   /*
   * Update ghost atom coordinates.
   *
   * Call on time steps for which no reneighboring is required. 
   */
   void Exchanger::updateGhosts()
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
               bufferPtr_->beginSendBlock(Buffer::GHOST);
               size = sendArray_(i, j).size();
               for (k = 0; k < size; ++k) {
                  bufferPtr_->packGhost(sendArray_(i, j)[k]);
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
                  bufferPtr_->unpackGhost(*atomPtr);
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
