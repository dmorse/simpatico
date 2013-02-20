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
#include <util/format/Dbl.h>
#include <util/global.h>

#include <algorithm>
#include <string>

#ifdef UTIL_DEBUG
//#define DDMD_EXCHANGER_DEBUG
#endif

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
      if (UTIL_ORTHOGONAL && !atomStoragePtr_->isCartesian()) {
         UTIL_THROW("Error: Coordinates not Cartesian on entry to exchange");
      } 
      if (!UTIL_ORTHOGONAL && atomStoragePtr_->isCartesian()) {
         UTIL_THROW("Error: Coordinates are Cartesian on entry to exchange");
      } 
      exchangeAtoms();
      exchangeGhosts();
   }

   /*
   * Initialize group ghost communication plan, before exchanging atoms.
   *
   * This method is called by exchangeAtoms, after computing plans for
   * exchanging atoms, based on their position, but before exchanging
   * atoms, and before clearing ghosts from previous exchange.
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
   * After calculating ghost communication plan for each group, clear 
   * the pointers to all ghost atoms in the group. The exchangeAtoms 
   * method will clear the actual ghost atoms from the AtomStorage.
   */
   template <int N>
   void Exchanger::initGroupGhostPlan (GroupStorage<N>& storage)
   {
      double coordinate;
      GroupIterator<N> groupIter;
      Atom* atomPtr;
      //int nEx[2];
      int nIn;
      int nOut;
      int i, j, k;
      bool choose;

      // Loop over groups
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
                  // nEx[j] = 0;
                  // Loop over atoms in group
                  for (k = 0; k < N; ++k) {
                     atomPtr = groupIter->atomPtr(k);
                     if (atomPtr) {
                        coordinate = atomPtr->position()[i];
                        if (atomPtr->isGhost()) {
                           if (j == 0) {
                              assert(inner_(i, j) > bound_(i, j));
                              if (coordinate < inner_(i, j)) {
                                 ++nOut;
                                 //if (coordinate < bound_(i, j)) {
                                 //   nEx[j] += 1;
                                 //}
                              }
                              if (coordinate > outer_(i, j)) {
                                 ++nIn;
                              }
                           } else { // if j = 1
                              assert(inner_(i, j) < bound_(i, j));
                              if (coordinate > inner_(i, j)) {
                                 ++nOut;
                                 //if (coordinate > bound_(i, j)) {
                                 //   nEx[j] += 1;
                                 //}
                              }
                              if (coordinate < outer_(i, j)) {
                                 ++nIn;
                              }
                           }
                        } else { // if atomPtr points to local atom
                           if (atomPtr->plan().exchange(i, j)) {
                              ++nOut;
                              //nEx[j] += 1;
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

               #if 0
               if (nEx[0] > 0 && nEx[1] > 0) {
                  UTIL_THROW("Group spanning 3 nodes");  
               }
               #endif

            } // end if multiProcessorDirection
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
   * Pack groups that contain atoms marked for exchange in this 
   * direction (direction i, j).
   *
   * Algorithm: Loop over groups. If the group contains one or
   * more atoms that are marked for exchange in direction i, j, 
   * pack the group for sending along. If the group is empty,
   * add it to emptyGroups, to mark it for later removal.
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
   * Remove empty groups from GroupStorage<N>.
   */
   template <int N>
   void Exchanger::removeEmptyGroups(GroupStorage<N>& storage,
                                     APArray< Group<N> >& emptyGroups)
   {
      int nEmpty = emptyGroups.size();
      #ifdef UTIL_DEBUG
      #ifdef DDMD_EXCHANGER_DEBUG
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

   /*
   * Unpack groups into bondStorage.
   */
   template <int N>
   void Exchanger::unpackGroups(GroupStorage<N>& storage)
   {
      Group<N>* newGroupPtr;
      Group<N>* oldGroupPtr;
      int groupId;

      bufferPtr_->beginRecvBlock();
      while (bufferPtr_->recvSize() > 0) {
         newGroupPtr = storage.newPtr();
         bufferPtr_->unpackGroup<N>(*newGroupPtr);
         groupId = newGroupPtr->id();
         oldGroupPtr = storage.find(groupId);
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
   #endif // endif ifdef UTIL_MPI

   /*
   * Set ghost communication flags for all atoms in incomplete groups.
   *
   * Precondition: This is called by exchangeAtoms after exchanging atoms 
   * and groups between neighboring processors. All ghosts are cleared.
   *
   * Algorithm: Loop over all Group<N> objects in the group storage. 
   * For each group, check if the group is incomplete, implying that one or
   * more atoms in the group are owned by another processor. If the group 
   * is incomplete, loop over 6 transfer directions. For each direction,
   * if the group is marked for sending in that direction, set the ghost
   * ghost communication flag for transfer in that direction for every 
   * local atom in the group. 
   *
   * Note: The algorithm assumes that the ghost communication flag for each
   * atom within a group whose atoms are divided among two or more processors 
   * will be set on the processor that owns the atom.
   */
   template <int N>
   void Exchanger::finishGroupGhostPlan(GroupStorage<N>& storage)
   {
      GroupIterator<N> groupIter;
      Atom* atomPtr;
      Plan* planPtr;
      int i, j, k, nAtom;

      // Loop over groups
      storage.begin(groupIter);
      for ( ; groupIter.notEnd(); ++groupIter) {

         #ifdef UTIL_DEBUG
         #ifdef DDMD_EXCHANGER_DEBUG
         // Validate group
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
            UTIL_THROW("Empty group");
         }
         #endif // ifdef DDMD_EXCHANGER_DEBUG
         #endif // ifdef UTIL_DEBUG

         // If this group is incomplete, set ghost flags for atoms 
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
                              planPtr = &atomPtr->plan();
                              if (!planPtr->ghost(i, j)) { 
                                 planPtr->setGhost(i, j);
                                 sendArray_(i,j).append(*atomPtr);
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
   * Exchange ownership of local atoms and groups, set ghost plan.
   *
   * Atomic coordinates must be in generalized coordinates on entry and 
   * exit if not UTIL_ORTHOGONAL, and Cartesian if UTIL_ORTHOGONAL.
   *
   * Algorithm:
   *
   *    - Loop over local atoms, set exchange and ghost communication
   *      flags for those beyond or near processor domain boundaries. 
   *
   *    - Add local atoms that will be retained by this processor but 
   *      sent as ghosts to appropriate send arrays.
   *
   *    - Call initGroupPlan each type of group (bond, angle, dihedral).
   *      initGroupPlan<N>{
   *
   *         For each group{
   *            - Set ghost communication flags for groups that span 
   *              or may span boundaries.
   *            - Clear pointers to ghost atoms in the group.
   *         }
   *      }
   *
   *      Clear all ghosts from the AtomStorage
   *
   *      For each transfer directions (i and j) {
   *  
   *         For each local atom {
   *            if marked for exchange(i, j) {
   *               if gridDimension[i] > 1 {
   *                  - add to sendAtoms array for removal
   *                  - pack into send buffer
   *               } else {
   *                  shift position to apply periodic b.c.
   *               }
   *            }
   *         }
   *
   *         if gridDimension[i] > 1 {
   *            - Call packGroups for each group type. 
   *              This packs groups containing atoms that are sent.
   *            - Remove exchanged atoms and empty groups
   *            - Send and receive data buffers
   *            for each atom in the receive buffer {
   *               - Unpack atom into AtomStorage
   *               - shift periodic boundary conditions
   *               - Determine if this is new home (or way station)
   *               - If atom is home, add to appropriate ghost arrays.
   *            }
   *            - Call unpackGroups for each group type.
   *         }
   * 
   *      } 
   *
   *    - Call finishGroupPlan<N> each type of group (bond, angle, dihedral).
   *      finishGroupPlan<N> {
   *         for each group{
   *            if group is incomplete{
   *               for each direction (i and j) {
   *                  if group is marked for ghost communication {
   *                     mark local atoms in group for ghost communication
   *                  }
   *               }
   *            }
   *         }
   *      }
   *
   *   }
   *
   *  Postconditions. Upon return:
   *     Each processor owns all atoms in its domain.
   *     Each processor owns all atoms containing one or more local atoms.
   *     Ghost plans are set for all local atoms.
   *     Send arrays contain local atoms marked for sending as ghosts.
   *     The AtomStorage contains no ghost atoms.
   */
   void Exchanger::exchangeAtoms()
   {
      stamp(START);
      Vector lengths = boundaryPtr_->lengths();
      double bound, inner, slabWidth;
      double coordinate, rshift;
      AtomIterator atomIter;
      Atom* atomPtr;
      Plan* planPtr;
      int i, j, jc, ip, jp, k, source, dest, nSend;
      int shift;
      bool isHome;
      bool isGhost;

      // Set domain and slab boundaries
      for (i = 0; i < Dimension; ++i) {
         if (UTIL_ORTHOGONAL) {
            slabWidth = pairCutoff_;
         } else {
            slabWidth = pairCutoff_/lengths[i];
         }
         for (j = 0; j < 2; ++j) {
            // j = 0 sends to lower coordinate i, bound is minimum
            // j = 1 sends to higher coordinate i, bound is maximum
            bound = domainPtr_->domainBound(i, j);
            bound_(i, j) = bound;
            if (j == 0) { // Communicate with lower index
               inner_(i,j) = bound + slabWidth;
               outer_(i, j)= bound - slabWidth;
            } else { // j == 1, communicate with upper index
               inner_(i, j) = bound - slabWidth;
               outer_(i, j) = bound + slabWidth;
            }
            sendArray_(i, j).clear();
         }
         if (domainPtr_->grid().dimension(i) > 1) {
            multiProcessorDirection_[i] = 1;
         } else {
            multiProcessorDirection_[i] = 0;
         }
      }

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
                     if (multiProcessorDirection_[i]) {
                        isHome = false;
                     }
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
                     if (multiProcessorDirection_[i]) {
                        isHome = false;
                     }
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
      #ifdef DDMD_EXCHANGER_DEBUG
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
            inner = inner_(i, jc);                 // inner bound upon receipt 
            shift = domainPtr_->shift(i, j);       // shift for periodic b.c.
            if (UTIL_ORTHOGONAL) {
               rshift = lengths[i]*double(shift);
            } else {
               rshift = double(shift);
            }

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
               coordinate = atomIter->position()[i];
               #ifdef DDMD_EXCHANGER_DEBUG
               {
                  bool choose;
                  if (j == 0) {
                     choose = (coordinate < bound);
                  } else {
                     choose = (coordinate > bound);
                  }
                  assert(choose == atomIter->plan().exchange(i, j));
               }
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

                     #ifdef UTIL_DEBUG
                     #ifdef DDMD_EXCHANGER_DEBUG 
                     assert(shift);
                     assert(coordinate > -1.0*fabs(rshift));
                     assert(coordinate <  2.0*fabs(rshift));
                     #endif
                     #endif

                     // Shift position if required by periodic b.c.
                     if (shift) {
                        atomIter->position()[i] += rshift;
                     }

                     #ifdef UTIL_DEBUG
                     coordinate = atomIter->position()[i];
                     assert(coordinate >= domainPtr_->domainBound(i, 0));
                     assert(coordinate < domainPtr_->domainBound(i, 1));
                     #endif

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

               // Pack groups that contain postmarked atoms.
               packGroups<2>(i, j, *bondStoragePtr_, emptyBonds_);
               #ifdef INTER_ANGLE
               packGroups<3>(i, j, *angleStoragePtr_, emptyAngles_);
               #endif
               #ifdef INTER_DIHEDRAL
               packGroups<4>(i, j, *dihedralStoragePtr_, emptyDihedrals_);
               #endif
               stamp(PACK_GROUPS);

               /*
               * Note: Removal cannot be done within above loops over atoms
               * and groups because element removal invalidates the iterators.
               */

               // Remove chosen atoms (listed in recvArray) from atomStorage
               nSend = sentAtoms_.size();
               for (k = 0; k < nSend; ++k) {
                  atomStoragePtr_->removeAtom(&sentAtoms_[k]);
               }
               stamp(REMOVE_ATOMS);
     
               // Remove empty groups
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
                  planPtr = &atomPtr->plan();
                  bufferPtr_->unpackAtom(*atomPtr);
                  atomStoragePtr_->addNewAtom();

                  if (shift) {
                     atomPtr->position()[i] += rshift;
                  }

                  #ifdef UTIL_DEBUG
                  // Check bounds on atom coordinate
                  coordinate = atomPtr->position()[i];
                  assert(coordinate > domainPtr_->domainBound(i, 0));
                  assert(coordinate < domainPtr_->domainBound(i, 1));

                  // Check ghost plan
                  assert(!planPtr->ghost(i, j));
                  if (j == 0) {
                     assert(planPtr->ghost(i, 1) == (coordinate > inner) );
                  } else {
                     assert(planPtr->ghost(i, 0) == (coordinate < inner) );
                  }
                  #endif

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

                  // If atom will stay, add to sendArrays for ghosts
                  if (isHome) {
                     for (ip = 0; ip < Dimension; ++ip) {
                        for (jp = 0; jp < 2; ++jp) {
                           if (planPtr->ghost(ip, jp)) {
                              sendArray_(ip, jp).append(*atomPtr);
                           }
                        }
                     }
                  }

               }
               assert(bufferPtr_->recvSize() == 0);
               stamp(UNPACK_ATOMS);

               // Unpack groups
               unpackGroups<2>(*bondStoragePtr_);
               #ifdef INTER_ANGLE
               unpackGroups<3>(*angleStoragePtr_);
               #endif
               #ifdef INTER_DIHEDRAL
               unpackGroups<4>(*dihedralStoragePtr_);
               #endif
               stamp(UNPACK_GROUPS);

            } // end if gridDimension > 1
            #endif // ifdef UTIL_MPI

         } // end for j (direction 0, 1)
      } // end for i (Cartesian index)

      /*
      * At this point:
      *    No ghost atoms exist.
      *    All atoms are on correct processor.
      *    No Groups are empty.
      *    All pointer to local atoms in Groups are set.
      *    All pointers to ghost atoms in Groups are null.
      */

      #ifdef UTIL_DEBUG
      #ifdef DDMD_EXCHANGER_DEBUG
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
      #endif // ifdef DDMD_EXCHANGER_DEBUG
      #endif // ifdef UTIL_DEBUG

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
      double  rshift;
      Atom* atomPtr;
      Atom* sendPtr;
      int i, j, jc, ip, jp, k, source, dest, shift, size;

      // Clear all receive arrays
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {
            recvArray_(i, j).clear();
         }
      }

      #ifdef UTIL_DEBUG
      double coordinate;
      #ifdef DDMD_EXCHANGER_DEBUG
      // Check send arrays
      {
         // Count local atoms marked for sending as ghosts.
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
         // Check consistency of counts with sendArray sizes.
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j < 2; ++j) {
               if (sendCounter(i, j) != sendArray_(i, j).size()) {
                  UTIL_THROW("Incorrect sendArray size");
               }
            }
         }
      }
      #endif // ifdef DDMD_EXCHANGER_DEBUG
      #endif // ifdef UTIL_DEBUG
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

            // Shift on receiving node for periodic b.c.s
            shift = domainPtr_->shift(i, j);
            if (UTIL_ORTHOGONAL) {
               rshift = lengths[i]*shift;
            } else {
               rshift = 1.0*shift;
            }

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
                     atomPtr->position()[i] += rshift;
                  }
                  atomStoragePtr_->addNewGhost();

                  #ifdef UTIL_DEBUG
                  // Validate shifted ghost coordinate 
                  coordinate = atomPtr->position()[i];
                  if (j == 0) {
                     assert(coordinate > bound_(i, 1));
                  } else {
                     assert(coordinate < bound_(i, 0));
                  }
                  #endif

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
                     atomPtr->position()[i] += rshift;
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
                  coordinate = atomPtr->position()[i];
                  if (j == 0) {
                     assert(coordinate > bound_(i, 1));
                  } else {
                     assert(coordinate < bound_(i, 0));
                  }
                  #endif

               }
               stamp(UNPACK_GHOSTS);

            }
            #endif // ifdef UTIL_MPI

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
      #ifdef DDMD_EXCHANGER_DEBUG
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
      if (!atomStoragePtr_->isCartesian()) {
         UTIL_THROW("Error: Coordinates not Cartesian on entry to update");
      } 

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
                     boundaryPtr_->applyShift(atomPtr->position(), i, shift);
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
                     boundaryPtr_->applyShift(atomPtr->position(), i, shift);
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
      int    i, j, k, source, dest, size;

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

   /*
   * Output statistics.
   */
   void Exchanger::outputStatistics(std::ostream& out, double time, int nStep)
   {
      // Precondition
      if (!domainPtr_->isMaster()) {
         UTIL_THROW("May be called only on domain master");
      }

      int nAtomTot = atomStoragePtr_->nAtomTotal();
      int nProc = 1;
      #ifdef UTIL_MPI
      nProc = domainPtr_->communicator().Get_size();
      #endif

      double ratio = double(nProc)/double(nStep*nAtomTot);

      out << std::endl;
      double AtomPlanT =  timer_.time(Exchanger::ATOM_PLAN);
      out << "AtomPlan              " << Dbl(AtomPlanT*ratio, 12, 6)
          << " sec   " << Dbl(AtomPlanT/time, 12, 6, true) << std::endl;
      double InitGroupPlanT =  timer_.time(Exchanger::INIT_GROUP_PLAN);
      out << "InitGroupPlan         " << Dbl(InitGroupPlanT*ratio, 12, 6)
          << " sec   " << Dbl(InitGroupPlanT/time, 12, 6, true) << std::endl;
      double ClearGhostsT =  timer_.time(Exchanger::CLEAR_GHOSTS);
      out << "ClearGhosts           " << Dbl(ClearGhostsT*ratio, 12, 6)
          << " sec   " << Dbl(ClearGhostsT/time, 12, 6, true) << std::endl;
      double PackAtomsT =  timer_.time(Exchanger::PACK_ATOMS);
      out << "PackAtoms             " << Dbl(PackAtomsT*ratio, 12, 6)
          << " sec   " << Dbl(PackAtomsT/time, 12, 6, true) << std::endl;
      double PackGroupsT =  timer_.time(Exchanger::PACK_GROUPS);
      out << "PackGroups            " << Dbl(PackGroupsT*ratio, 12, 6)
          << " sec   " << Dbl(PackGroupsT/time, 12, 6, true) << std::endl;
      double RemoveAtomsT =  timer_.time(Exchanger::REMOVE_ATOMS);
      out << "RemoveAtoms           " << Dbl(RemoveAtomsT*ratio, 12, 6)
          << " sec   " << Dbl(RemoveAtomsT/time, 12, 6, true) << std::endl;
      double RemoveGroupsT =  timer_.time(Exchanger::REMOVE_GROUPS);
      out << "RemoveGroups          " << Dbl(RemoveGroupsT*ratio, 12, 6)
          << " sec   " << Dbl(RemoveGroupsT/time, 12, 6, true) << std::endl;
      double SendRecvAtomsT =  timer_.time(Exchanger::SEND_RECV_ATOMS);
      out << "SendRecvAtoms         " << Dbl(SendRecvAtomsT*ratio, 12, 6)
          << " sec   " << Dbl(SendRecvAtomsT/time, 12, 6, true) << std::endl;
      double UnpackAtomsT =  timer_.time(Exchanger::UNPACK_ATOMS);
      out << "UnpackAtoms           " << Dbl(UnpackAtomsT*ratio, 12, 6)
          << " sec   " << Dbl(UnpackAtomsT/time, 12, 6, true) << std::endl;
      double UnpackGroupsT =  timer_.time(Exchanger::UNPACK_GROUPS);
      out << "UnpackGroups          " << Dbl(UnpackGroupsT*ratio, 12, 6)
          << " sec   " << Dbl(UnpackGroupsT/time, 12, 6, true) << std::endl;
      double FinishGroupPlanT =  timer_.time(Exchanger::FINISH_GROUP_PLAN);
      out << "FinishGroupPlan       " << Dbl(FinishGroupPlanT*ratio, 12, 6)
          << " sec   " << Dbl(FinishGroupPlanT/time, 12, 6, true) << std::endl;
      double SendArraysT =  timer_.time(Exchanger::INIT_SEND_ARRAYS);
      out << "SendArrays            " << Dbl(SendArraysT*ratio, 12, 6)
          << " sec   " << Dbl(SendArraysT/time, 12, 6, true) << std::endl;
      double PackGhostsT =  timer_.time(Exchanger::PACK_GHOSTS);
      out << "PackGhosts            " << Dbl(PackGhostsT*ratio, 12, 6)
          << " sec   " << Dbl(PackGhostsT/time, 12, 6, true) << std::endl;
      double SendRecvGhostsT =  timer_.time(Exchanger::SEND_RECV_GHOSTS);
      out << "SendRecvGhosts        " << Dbl(SendRecvGhostsT*ratio, 12, 6)
          << " sec   " << Dbl(SendRecvGhostsT/time, 12, 6, true) << std::endl;
      double UnpackGhostsT =  timer_.time(Exchanger::UNPACK_GHOSTS);
      out << "UnpackGhosts          " << Dbl(UnpackGhostsT*ratio, 12, 6)
          << " sec   " << Dbl(UnpackGhostsT/time, 12, 6, true) << std::endl;
      double FindGroupGhostsT =  timer_.time(Exchanger::FIND_GROUP_GHOSTS);
      out << "FindGroupGhosts       " << Dbl(FindGroupGhostsT*ratio, 12, 6)
          << " sec   " << Dbl(FindGroupGhostsT/time, 12, 6, true) << std::endl;
      double PackUpdateT =  timer_.time(Exchanger::PACK_UPDATE);
      out << "PackUpdate            " << Dbl(PackUpdateT*ratio, 12, 6)
          << " sec   " << Dbl(PackUpdateT/time, 12, 6, true) << std::endl;
      double SendRecvUpdateT =  timer_.time(Exchanger::SEND_RECV_UPDATE);
      out << "SendRecvUpdate        " << Dbl(SendRecvUpdateT*ratio, 12, 6)
          << " sec   " << Dbl(SendRecvUpdateT/time, 12, 6, true) << std::endl;
      double UnpackUpdateT =  timer_.time(Exchanger::UNPACK_UPDATE);
      out << "UnpackUpdate          " << Dbl(UnpackUpdateT*ratio, 12, 6)
          << " sec   " << Dbl(UnpackUpdateT/time, 12, 6, true) << std::endl;
      double LocalUpdateT =  timer_.time(Exchanger::LOCAL_UPDATE);
      out << "LocalUpdate           " << Dbl(LocalUpdateT*ratio, 12, 6)
          << " sec   " << Dbl(LocalUpdateT/time, 12, 6, true) << std::endl;
      out << std::endl;

   }
   #endif // ifdef UTIL_MPI

}
#endif
