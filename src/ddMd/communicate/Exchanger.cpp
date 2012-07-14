#ifndef DDMD_EXCHANGER_CPP
#define DDMD_EXCHANGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
#define EXCHANGER_GHOST

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   Exchanger::Exchanger()
    : sendArray_(),
      recvArray_(),
      #ifdef UTIL_MPI
      emptyBonds_(),
      #ifdef INTER_ANGLE
      emptyAngles_(),
      #endif
      #ifdef INTER_DIHEDRAL
      emptyDihedrals_(),
      #endif
      #endif
      bound_(),
      inner_(),
      outer_(),
      multiProcessorDirection_(),
      boundaryPtr_(0),
      domainPtr_(0),
      atomStoragePtr_(0),
      bondStoragePtr_(0),
      #ifdef INTER_ANGLE 
      angleStoragePtr_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_(0),
      #endif
      bufferPtr_(0),
      pairCutoff_(-1.0),
      timer_(Exchanger::NTime)
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

      #ifdef UTIL_MPI
      emptyBonds_.reserve(sendRecvCapacity);
      #ifdef INTER_ANGLE
      emptyAngles_.reserve(sendRecvCapacity);
      #endif
      #ifdef INTER_DIHEDRAL
      emptyDihedrals_.reserve(sendRecvCapacity);
      #endif
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
            if (multiProcessorDirection_[i]) {
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

   #ifdef UTIL_MPI
   /*
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
      int nEmpty = emptyGroups.size();
      #ifdef UTIL_DEBUG
      #ifdef EXCHANGER_DEBUG
      // Confirm that groups are actually empty
      Atom* atomPtr;
      int   atomId;
      for (int k = 0; k < nEmpty; ++k) {
         for (int m = 0; m < N; ++m) {
            atomId = emptyGroups[k].atomId(m);
            atomPtr = emptyGroups[k].atomPtr(m);
            assert(atomPtr == 0);
            assert(atomStoragePtr_->find(atomId) == 0);
         }
      }
      #endif
      #endif
      for (int k = 0; k < nEmpty; ++k) {
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
   #endif

   template <int N>
   void Exchanger::finishGroupGhostPlan(GroupStorage<N>& storage)
   {
      GroupIterator<N> groupIter;
      Atom* atomPtr;
      Plan* planPtr;
      int i, j, k, nAtom;

      // Set ghost communication flags for atoms in incomplete bonds
      storage.begin(groupIter);
      for ( ; groupIter.notEnd(); ++groupIter) {

         #ifdef UTIL_DEBUG
         #ifdef EXCHANGER_DEBUG
         int atomId;
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
               if (multiProcessorDirection_[i]) {
                  for (j = 0; j < 2; ++j) {
                     if (groupIter->plan().ghost(i, j)) {
                        for (k = 0; k < N; ++k) {
                           atomPtr = groupIter->atomPtr(k);
                           if (atomPtr) {
                              assert(!atomPtr->isGhost());
                              #ifdef EXCHANGER_GHOST
                              planPtr = &atomPtr->plan();
                              if (!planPtr->ghost(i, j)) { 
                                 planPtr->setGhost(i, j);
                                 sendArray_(i,j).append(*atomPtr);
                              }
                              #else
                              atomPtr->plan().setGhost(i, j);
                              #endif
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
      stamp(START);
      Vector lengths = boundaryPtr_->lengths();
      double bound, inner;
      double coordinate;
      AtomIterator atomIter;
      Atom* atomPtr;
      Plan* planPtr;
      int i, j, jc, ip, jp, k, source, dest, nSend;
      int shift;
      bool isHome;
      bool isGhost;

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
            sendArray_(i, j).clear();
         }
         if (domainPtr_->grid().dimension(i) > 1) {
            multiProcessorDirection_[i] = 1;
         } else {
            multiProcessorDirection_[i] = 0;
         }
      }

      /**
      * In this function:
      *  - sendArray_ is filled with ptrs to ghosts for transfer
      *  - recvArray_ is filled with ptrs to local atoms for removal
      * In exchangeGhosts(), recvArray is cleared, and then refilled 
      * with pointers to ghosts as they are received (hence the name).
      */

      // Compute communication plan for every local atom
      atomStoragePtr_->begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {

         planPtr = &atomIter->plan();
         planPtr->clearFlags();
         isHome  = true;
         isGhost = false;

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
                     planPtr->setExchange(i, j);
                     #ifdef EXCHANGER_GHOST
                     if (multiProcessorDirection_[i]) {
                        isHome = false;
                     }
                     #endif
                     if (coordinate > outer_(i, j)) {
                        planPtr->setGhost(i, jc);
                        isGhost = true;
                     }
                  } else {
                     if (coordinate < inner_(i, j)) {
                        planPtr->setGhost(i, j);
                        isGhost = true;
                     }
                  }
               } else { // j == 1, communicate with upper index
                  if (coordinate > bound_(i, j)) {
                     planPtr->setExchange(i, j);
                     #ifdef EXCHANGER_GHOST
                     if (multiProcessorDirection_[i]) {
                        isHome = false;
                     }
                     #endif
                     if (coordinate < outer_(i, j)) {
                        planPtr->setGhost(i, jc);
                        isGhost = true;
                     }
                  } else {
                     if (coordinate > inner_(i, j)) {
                        planPtr->setGhost(i, j);
                        isGhost = true;
                     }
                  }
               }
   
            } // end for j
         } // end for i

         #ifdef EXCHANGER_GHOST
         // Add atoms that will be retained by this processor,
         // but will be communicated as ghosts to sendArray_
         if (isGhost && isHome) {
            for (i = 0; i < Dimension; ++i) {
               for (j = 0; j < 2; ++j) {
                  if (planPtr->ghost(i, j)) {
                     sendArray_(i, j).append(*atomIter);
                  }
               }
            }
         }
         #endif

      } // end atom loop, end compute plan
      stamp(ATOM_PLAN);

      /*
      * Find groups that span boundaries (uses information about ghosts).
      * Clear pointers to ghosts in each Group after inspecting the Group.
      *
      * Atoms in Group<N> objects that span boundaries, and are incomplete
      * after atom migration, are marked for sending as ghosts in the
      * finishGroupGhostPlan<N> function, further below.
      */
      initGroupGhostPlan<2>(*bondStoragePtr_);      // bonds
      #ifdef INTER_ANGLE
      initGroupGhostPlan<3>(*angleStoragePtr_);     // angles
      #endif
      #ifdef INTER_DIHEDRAL
      initGroupGhostPlan<4>(*dihedralStoragePtr_);  // dihedrals
      #endif
      stamp(INIT_GROUP_PLAN);

      // Clear all ghost atoms from AtomStorage
      atomStoragePtr_->clearGhosts();
      stamp(CLEAR_GHOSTS);

      #ifdef UTIL_DEBUG
      #ifdef EXCHANGER_DEBUG
      int nAtomTotal;
      atomStoragePtr_->computeNAtomTotal(domainPtr_->communicator());
      int myRank = domainPtr_->gridRank();
      if (myRank == 0) {
         nAtomTotal = atomStoragePtr_->nAtomTotal();
      }
      bondStoragePtr_->isValid(*atomStoragePtr_, 
                               domainPtr_->communicator(), false);
      #endif
      #endif

      // Cartesian directions for exchange (0=x, 1=y, 2=z)
      for (i = 0; i < Dimension; ++i) {

         // Transmission direction
         // j = 0 sends to processor with lower grid coordinate i
         // j = 1 sends to processor with higher grid coordinate i
         for (j = 0; j < 2; ++j) {

            // Index for conjugate (reverse) direction
            if (j == 0) jc = 1;
            if (j == 1) jc = 0;

            source = domainPtr_->sourceRank(i, j); // rank to receive from
            dest = domainPtr_->destRank(i, j);     // rank to send to
            bound = domainPtr_->domainBound(i, j); // bound for send
            inner = inner_(i, jc);             // inner bound upon receipt 
            shift = domainPtr_->shift(i, j);   // shift for periodic b.c.

            #ifdef UTIL_MPI
            if (multiProcessorDirection_[i]) {
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::ATOM);
            }
            #endif

            // Choose atoms for sending, pack and mark for removal.
            sentAtoms_.clear();
            atomStoragePtr_->begin(atomIter);
            for ( ; atomIter.notEnd(); ++atomIter) {

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

               if (atomIter->plan().exchange(i, j)) {

                  #ifdef UTIL_MPI
                  if (multiProcessorDirection_[i]) {

                     sentAtoms_.append(*atomIter);
                     bufferPtr_->packAtom(*atomIter);

                  } else 
                  #endif
                  {

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
            stamp(PACK_ATOMS);

            #ifdef UTIL_MPI
            // Send and receive only if processor grid dimension(i) > 1
            if (multiProcessorDirection_[i]) {

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
               stamp(PACK_GROUPS);

               /*
               * Note: Removal cannot be done within the above loops over 
               * atoms and groups because element removal invalidates a
               * PArray iterator.
               */

               // Remove chosen atoms (listed in recvArray) from atomStorage
               nSend = sentAtoms_.size();
               for (k = 0; k < nSend; ++k) {
                  atomStoragePtr_->removeAtom(&sentAtoms_[k]);
               }
               stamp(REMOVE_ATOMS);
     
               // Remove empty bonds from bondStorage.
               removeEmptyGroups<2>(*bondStoragePtr_, emptyBonds_);
               #ifdef INTER_ANGLE
               removeEmptyGroups<3>(*angleStoragePtr_, emptyAngles_);
               #endif
               #ifdef INTER_DIHEDRAL
               removeEmptyGroups<4>(*dihedralStoragePtr_, emptyDihedrals_);
               #endif
               stamp(REMOVE_GROUPS);

               // Send to processor dest and receive from processor source
               bufferPtr_->sendRecv(domainPtr_->communicator(), 
                                    source, dest);
               stamp(SEND_RECV_ATOMS);

               // Unpack atoms into atomStorage
               bufferPtr_->beginRecvBlock();
               while (bufferPtr_->recvSize() > 0) {

                  atomPtr = atomStoragePtr_->newAtomPtr();
                  bufferPtr_->unpackAtom(*atomPtr);
                  atomStoragePtr_->addNewAtom();

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

                  #ifdef EXCHANGER_GHOST
                  planPtr = &atomPtr->plan();

                  // Determine if new atom will stay on this processor.
                  isHome = true;
                  if (i < Dimension - 1) {
                     for (ip = i + 1; ip < Dimension; ++ip) {
                        if (multiProcessorDirection_[ip]) {
                           for (jp = 0; jp < 2; ++jp) {
                              if (planPtr->exchange(ip, jp)) {
                                 isHome = false;
                              }
                           }
                        }
                    }
                  }

                  // If atom will stay, add to ghost sendArrays
                  if (isHome) {
                     for (ip = 0; ip < Dimension; ++ip) {
                        for (jp = 0; jp < 2; ++jp) {
                           if (planPtr->ghost(ip, jp)) {
                              sendArray_(ip, jp).append(*atomPtr);
                           }
                        }
                     }
                  }
                  #endif

               }
               assert(bufferPtr_->recvSize() == 0);
               stamp(UNPACK_ATOMS);

               // Unpack bonds into bondStorage.
               unpackGroups<2>(*bondStoragePtr_);
               #ifdef INTER_ANGLE
               unpackGroups<3>(*angleStoragePtr_);
               #endif
               #ifdef INTER_DIHEDRAL
               unpackGroups<4>(*dihedralStoragePtr_);
               #endif
               stamp(UNPACK_GROUPS);

            } // end if gridDimension > 1
            #endif

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
      // Validity checks
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

      // Set ghost communication flags for atoms in incomplete groups
      finishGroupGhostPlan<2>(*bondStoragePtr_);
      #ifdef INTER_ANGLE
      finishGroupGhostPlan<3>(*angleStoragePtr_);
      #endif
      #ifdef INTER_DIHEDRAL
      finishGroupGhostPlan<4>(*dihedralStoragePtr_);
      #endif

      stamp(FINISH_GROUP_PLAN);
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
      for ( ; groupIter.notEnd(); ++groupIter) {
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
      stamp(START);

      // Preconditions
      assert(bufferPtr_->isInitialized());
      assert(domainPtr_->isInitialized());
      assert(domainPtr_->hasBoundary());
      if (atomStoragePtr_->nGhost() != 0) {
         UTIL_THROW("atomStoragePtr_->nGhost() != 0");
      }

      Vector  lengths = boundaryPtr_->lengths();
      double  bound, inner;
      Atom* atomPtr;
      Atom* sendPtr;
      int i, j, jc, ip, jp, k, source, dest, shift, size;

      // Clear all send and receive arrays
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {
            recvArray_(i, j).clear();
         }
      }

      #ifndef EXCHANGER_GHOST
      AtomIterator  localIter;
      atomStoragePtr_->begin(localIter);
      for ( ; localIter.notEnd(); ++localIter) {
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j < 2; ++j) {
               if (localIter->plan().ghost(i, j)) {
                   sendArray_(i, j).append(*localIter);
               }
            }
         }
      }
      #else
      #ifdef UTIL_DEBUG
      #ifdef EXCHANGER_DEBUG
      // Add local atoms to all appropriate send arrays
      {
         FMatrix<int, Dimension, 2> sendCounter;
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j < 2; ++j) {
               sendCounter(i, j) = 0;
            }
         }
         AtomIterator  localIter;
         Plan* planPtr;
         atomStoragePtr_->begin(localIter);
         for ( ; localIter.notEnd(); ++localIter) {
            planPtr = &localIter->plan();
            for (i = 0; i < Dimension; ++i) {
               for (j = 0; j < 2; ++j) {
                  if (planPtr->ghost(i, j)) {
                      ++sendCounter(i, j);
                  }
               }
            }
         }
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j < 2; ++j) {
               if (sendCounter(i, j) != sendArray_(i, j).size()) {
                  UTIL_THROW("Incorrect sendArray size");
               }
            }
         }
      }
      #endif
      #endif
      #endif
      stamp(INIT_SEND_ARRAYS);

      // Cartesian directions
      for (i = 0; i < Dimension; ++i) {

         // Transmit directions
         for (j = 0; j < 2; ++j) {

            // j = 0: Send ghosts near minimum bound to lower coordinate
            // j = 1: Send ghosts near maximum bound to higher coordinate

            // Set index for reverse direction
            if (j == 0) jc = 1;
            if (j == 1) jc = 0;

            bound = bound_(i, j);
            inner = inner_(i, j);
       
            // Shift on receiving node for periodic b.c.s
            shift = domainPtr_->shift(i, j);

            #ifdef UTIL_MPI
            if (multiProcessorDirection_[i]) {
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::GHOST);
            }
            #endif

            // Pack atoms in sendArray_(i, j)
            size = sendArray_(i, j).size();
            for (k = 0; k < size; ++k) {

               sendPtr = &sendArray_(i, j)[k];

               #ifdef UTIL_MPI
               if (multiProcessorDirection_[i]) {

                  // If grid dimension > 1, pack atom for sending 
                  bufferPtr_->packGhost(*sendPtr);

               } else 
               #endif
               {  // if grid dimension == 1

                  // Make a ghost copy of local atom on this processor
                  atomPtr = atomStoragePtr_->newGhostPtr();
                  recvArray_(i, j).append(*atomPtr);
                  atomPtr->setId(sendPtr->id());
                  atomPtr->setTypeId(sendPtr->typeId());
                  atomPtr->plan().setFlags(sendPtr->plan().flags());
                  atomPtr->position() = sendPtr->position();
                  if (shift) {
                     atomPtr->position()[i] += shift * lengths[i];
                  }
                  atomStoragePtr_->addNewGhost();

                  #ifdef UTIL_DEBUG
                  // Validate shifted ghost coordinate 
                  if (j == 0) {
                     assert(atomPtr->position()[i] > bound_(i, 1));
                  } else {
                     assert(atomPtr->position()[i] < bound_(i, 0));
                  }
                  #endif

                  // Prevent ghost copy from being re-copied within
                  // following loop over ghosts on same processor.
                  // atomPtr->plan().clearGhost(i, j);
 
                  // Add to send arrays for any remaining directions
                  if (i < Dimension - 1) {
                     for (ip = i + 1; ip < Dimension; ++ip) {
                        for (jp = 0; jp < 2; ++jp) {
                           if (atomPtr->plan().ghost(ip, jp)) {
                              sendArray_(ip, jp).append(*atomPtr);
                           }
                        }
                     }
                  }

               }

            }
            stamp(PACK_GHOSTS);

            #ifdef UTIL_MPI
            // Send and receive buffers
            if (multiProcessorDirection_[i]) {

               bufferPtr_->endSendBlock();

               source = domainPtr_->sourceRank(i, j);
               dest   = domainPtr_->destRank(i, j);
               bufferPtr_->sendRecv(domainPtr_->communicator(), source, dest);
               stamp(SEND_RECV_GHOSTS);

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

                  // Add to send arrays for remaining directions
                  if (i < Dimension - 1) {
                     for (ip = i + 1; ip < Dimension; ++ip) {
                        for (jp = 0; jp < 2; ++jp) {
                           if (atomPtr->plan().ghost(ip, jp)) {
                              sendArray_(ip, jp).append(*atomPtr);
                           }
                        }
                     }
                  }

                  #ifdef UTIL_DEBUG
                  // Validate ghost coordinate on the receiving processor.
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
               stamp(UNPACK_GHOSTS);

            }
            #endif

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
      bondStoragePtr_->isValid(*atomStoragePtr_, 
                               domainPtr_->communicator(), true);
      #endif
      #endif

      stamp(FIND_GROUP_GHOSTS);
   }

   /*
   * Update ghost atom coordinates.
   *
   * Call on time steps for which no reneighboring is required. 
   */
   void Exchanger::update()
   {
      stamp(START);
      Vector lengths = boundaryPtr_->lengths();
      Atom*  atomPtr;
      int    i, j, k, source, dest, size, shift;

      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {

            // Shift on receiving processor for periodic boundary conditions
            shift = domainPtr_->shift(i, j);

            if (multiProcessorDirection_[i]) {

               // Pack ghost positions for sending
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::UPDATE);
               size = sendArray_(i, j).size();
               for (k = 0; k < size; ++k) {
                  bufferPtr_->packUpdate(sendArray_(i, j)[k]);
               }
               bufferPtr_->endSendBlock();
               stamp(PACK_UPDATE);
  
               // Send and receive buffers
               source = domainPtr_->sourceRank(i, j);
               dest   = domainPtr_->destRank(i, j);
               bufferPtr_->sendRecv(domainPtr_->communicator(), source, dest);
               stamp(SEND_RECV_UPDATE);
   
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
               stamp(UNPACK_UPDATE);

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
               stamp(LOCAL_UPDATE);

            }

         } // transmit direction j = 0, 1

      } // Cartesian direction i = 0, ..., Dimension - 1

   }

   /*
   * Update ghost atom forces.
   *
   * Call on time steps for which no reneighboring is required,
   * if reverse communication is enabled.
   */
   void Exchanger::reverseUpdate()
   {
      stamp(START);
      Vector lengths = boundaryPtr_->lengths();
      Atom*  atomPtr;
      int    i, j, k, source, dest, size, shift;

      for (i = Dimension - 1; i >= 0; --i) {
         for (j = 1; j >= 0; --j) {

            if (multiProcessorDirection_[i]) {

               // Pack ghost forces for sending
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::FORCE);
               size = recvArray_(i, j).size();
               for (k = 0; k < size; ++k) {
                  bufferPtr_->packForce(recvArray_(i, j)[k]);
               }
               bufferPtr_->endSendBlock();
               stamp(PACK_FORCE);
  
               // Send and receive buffers (reverse direction)
               source  = domainPtr_->destRank(i, j);
               dest    = domainPtr_->sourceRank(i, j);
               bufferPtr_->sendRecv(domainPtr_->communicator(), 
                                    source, dest);
               stamp(SEND_RECV_FORCE);
   
               // Unpack ghost forces
               bufferPtr_->beginRecvBlock();
               size = sendArray_(i, j).size();
               for (k = 0; k < size; ++k) {
                  atomPtr = &sendArray_(i, j)[k];
                  bufferPtr_->unpackForce(*atomPtr);
               }
               stamp(UNPACK_FORCE);

            } else {

               // If grid().dimension(i) == 1, then copy forces of atoms
               // listed in sendArray to those listed in the recvArray.

               size = recvArray_(i, j).size();
               assert(size == sendArray_(i, j).size());
               for (k = 0; k < size; ++k) {
                  atomPtr = &sendArray_(i, j)[k];
                  atomPtr->force() += recvArray_(i, j)[k].force();
               }
               stamp(LOCAL_FORCE);

            }

         } // transmit direction j = 1 or 0 

      } // Cartesian direction i 

   }
   #endif

}
#endif
