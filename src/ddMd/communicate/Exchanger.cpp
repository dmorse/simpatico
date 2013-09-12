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
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Group.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/storage/GroupExchanger.h>
#include <ddMd/storage/GroupIterator.h>
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
      bound_(),
      inner_(),
      outer_(),
      gridFlags_(),
      boundaryPtr_(0),
      domainPtr_(0),
      atomStoragePtr_(0),
      groupExchangers_(),
      bufferPtr_(0),
      pairCutoff_(-1.0),
      timer_(Exchanger::NTime)
   {  groupExchangers_.reserve(8); }

   /*
   * Destructor.
   */
   Exchanger::~Exchanger()
   {}

   /*
   * Set pointers to associated objects.
   */
   void Exchanger::associate(const Domain& domain, const Boundary& boundary,
                             AtomStorage& atomStorage, Buffer& buffer)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;
      atomStoragePtr_ = &atomStorage;
      bufferPtr_ = &buffer;
   }

   /*
   * Add a GroupExchanger to an internal list.
   */
   void Exchanger::addGroupExchanger(GroupExchanger& groupExchanger)
   {  groupExchangers_.append(groupExchanger); }

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
      if (atomStoragePtr_->isCartesian()) {
         UTIL_THROW("Error: Coordinates are Cartesian on entry to exchange");
      }
      exchangeAtoms();
      exchangeGhosts();
   }

   /*
   * Exchange ownership of local atoms and groups, set ghost plan.
   *
   * Atomic coordinates must be in scaled / generalized form on entry.
   *
   * Algorithm:
   *
   *    - Loop over local atoms, set exchange and ghost communication
   *      flags for those beyond or near processor domain boundaries.
   *
   *    - Add local atoms that will be retained by this processor but
   *      sent as ghosts to appropriate send arrays.
   *
   *    - Call initGroupGhostPlan each type of group (bond, angle, dihedral).
   *      initGroupGhostPlan<N>{
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
         slabWidth = pairCutoff_/lengths[i];
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
            gridFlags_[i] = 1;
         } else {
            gridFlags_[i] = 0;
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
                     if (gridFlags_[i]) {
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
                     if (gridFlags_[i]) {
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
      * Set ghost communication flags for groups that span boundaries and
      * clear pointers to ghosts in each Group after inspecting the Group.
      * This function information uses old coordinate information about 
      * ghosts atoms, and so must be called before all ghosts are cleared.
      * 
      * Algorithm:: GroupExchanger::markSpanningGroups loops over all groups,
      * Because it is called before ghosts are cleared, all groups should 
      * be complete, but coordinate values for ghosts are values obtained 
      * from the previous update, and so are outdated by one step. (The only
      * exception is the first step of a simulation, which is discussed
      * separately below). For each group, it loops over all directions 
      * (i, j), with i=0,..,Dimension-1 and j=0, 1. For each direction it 
      * determines which atoms are or might be "inside" or "outside" of the 
      * boundary. This is determining by comparing atom coordinate i to 
      * inner_(i, j) and outer_(i, j), which are scaled coordinates somehwat 
      * "inside" and "outside" of the true domain boundary. If coordinate i
      * is "outside" of inner_(i, j), the atom is considered potentially 
      * inside the boundary. If a coordinate is "inside" of outer_(i, j), 
      * it is considered potentially outside the boundary. If a coordinate 
      * is between inner_(i, j) and outer_(i, j), it is considered both 
      * inside or outside, to reflect the uncertainty that arises from 
      * the use of outdated ghost positions. 
      *
      * If a group has atoms both "inside" and "outside" boundary (i, j),
      * the group is said to "span" the boundary, and ghost communication 
      * flag (i, j) is set for the Group. Otherwise, this ghost communication 
      * flag for the Group cleared. After finishing inspection of a Group, 
      * all pointers to ghost atoms within the Group are cleared.
      *
      * If a Group is incomplete when this function is called, Group ghost
      * flags are set for all directions. This should only happen on the
      * first time step of a simulation. 
      *
      * Why this necessary: After all local atoms have been exchanged, atoms
      * that belong to incomplete groups are marked for sending as ghosts 
      * in each direction for which the group spans a boundary. This is
      * done in the GroupExchanger::markGhosts function at the end of this
      * function.
      *
      * The function interfaces defined in GroupExchanger are implemented 
      * in the GroupStorage<int N> class template. 
      */

      // Set ghost communication flags for groups (see above)
      for (k = 0; k < groupExchangers_.size(); ++k) {
         groupExchangers_[k].markSpanningGroups(bound_, inner_, outer_,
                                                gridFlags_);
      }
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
      for (k = 0; k < groupExchangers_.size(); ++k) {
         groupExchangers_[k].isValid(*atomStoragePtr_,
                                     domainPtr_->communicator(), false);
      }
      #endif
      #endif

      // Cartesian directions for exchange (0=x, 1=y, 2=z)
      for (i = 0; i < Dimension; ++i) {

         // Transmission direction
         // j = 0 sends down to processor with lower grid coordinate i
         // j = 1 sends up   to processor with higher grid coordinate i
         for (j = 0; j < 2; ++j) {

            // Index for conjugate (reverse) direction
            if (j == 0) jc = 1;
            if (j == 1) jc = 0;

            source = domainPtr_->sourceRank(i, j); // rank to receive from
            dest = domainPtr_->destRank(i, j);     // rank to send to
            bound = domainPtr_->domainBound(i, j); // bound for send
            inner = inner_(i, jc);                 // inner bound upon receipt
            shift = domainPtr_->shift(i, j);       // shift for periodic b.c.
            rshift = double(shift);

            #ifdef UTIL_MPI
            if (gridFlags_[i]) {
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
                  if (gridFlags_[i]) {
                     sentAtoms_.append(*atomIter);
                     atomIter->packAtom(*bufferPtr_);
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

            /*
            * Notes:
            *
            * (1) Removal of atoms cannot be done within the atom packing
            * loop because element removal would invalidate the atom 
            * iterator.
            *
            * (2) Groups must be packed for sending before atoms are removed
            * because the algorithm for identifying which groups to send 
            * references pointers to associated atoms.
            */

            #ifdef UTIL_MPI
            // Send and receive buffers only if processor grid dimension(i) > 1
            if (gridFlags_[i]) {

               // End atom send block
               bufferPtr_->endSendBlock();

               // Pack groups that contain atoms marked for exchange.
               // Remove empty groups from each GroupStorage.
               for (k = 0; k < groupExchangers_.size(); ++k) {
                  groupExchangers_[k].pack(i, j, *bufferPtr_);
               }
               stamp(PACK_GROUPS);

               // Remove chosen atoms (from sentAtoms) from atomStorage
               nSend = sentAtoms_.size();
               for (k = 0; k < nSend; ++k) {
                  atomStoragePtr_->removeAtom(&sentAtoms_[k]);
               }
               stamp(REMOVE_ATOMS);

               // Send to processor dest and receive from processor source
               bufferPtr_->sendRecv(domainPtr_->communicator(),
                                    source, dest);
               stamp(SEND_RECV_ATOMS);

               // Unpack atoms into atomStorage
               bufferPtr_->beginRecvBlock();
               while (bufferPtr_->recvSize() > 0) {

                  atomPtr = atomStoragePtr_->newAtomPtr();
                  planPtr = &atomPtr->plan();
                  atomPtr->unpackAtom(*bufferPtr_);
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
                        if (gridFlags_[ip]) {
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
               bufferPtr_->endRecvBlock();
               assert(bufferPtr_->recvSize() == 0);
               stamp(UNPACK_ATOMS);

               // Unpack groups
               for (k = 0; k < groupExchangers_.size(); ++k) {
                  groupExchangers_[k].unpack(*bufferPtr_, *atomStoragePtr_);
               }
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
      for (k = 0; k < groupExchangers_.size(); ++k) {
         groupExchangers_[k].isValid(*atomStoragePtr_,
                                      domainPtr_->communicator(), false);
      }
      #endif // ifdef DDMD_EXCHANGER_DEBUG
      #endif // ifdef UTIL_DEBUG

      // Set ghost communication flags for atoms in incomplete groups
      for (k = 0; k < groupExchangers_.size(); ++k) {
         groupExchangers_[k].markGhosts(*atomStoragePtr_, sendArray_,
                                        gridFlags_);
      }
      stamp(FINISH_GROUP_PLAN);
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
            rshift = 1.0*shift;

            #ifdef UTIL_MPI
            if (gridFlags_[i]) {
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::GHOST);
            }
            #endif

            // Pack atoms in sendArray_(i, j)
            size = sendArray_(i, j).size();
            for (k = 0; k < size; ++k) {

               sendPtr = &sendArray_(i, j)[k];

               #ifdef UTIL_MPI
               if (gridFlags_[i]) {
                  // If grid dimension > 1, pack atom for sending
                  sendPtr->packGhost(*bufferPtr_);
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
            if (gridFlags_[i]) {

               bufferPtr_->endSendBlock();

               source = domainPtr_->sourceRank(i, j);
               dest   = domainPtr_->destRank(i, j);
               bufferPtr_->sendRecv(domainPtr_->communicator(), source, dest);
               stamp(SEND_RECV_GHOSTS);

               // Unpack ghosts and add to recvArray
               bufferPtr_->beginRecvBlock();
               while (bufferPtr_->recvSize() > 0) {

                  atomPtr = atomStoragePtr_->newGhostPtr();
                  atomPtr->unpackGhost(*bufferPtr_);
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
               bufferPtr_->endRecvBlock();
               stamp(UNPACK_GHOSTS);

            }
            #endif // ifdef UTIL_MPI

         } // end for transmit direction j = 0, 1

      } // end for Cartesian index i


      // Find ghost atoms for all incomplete groups
      for (k = 0; k < groupExchangers_.size(); ++k) {
         groupExchangers_[k].findGhosts(*atomStoragePtr_);
      }

      #ifdef UTIL_DEBUG
      #ifdef DDMD_EXCHANGER_DEBUG
      atomStoragePtr_->isValid();
      for (k = 0; k < groupExchangers_.size(); ++k) {
         groupExchangers_[k].isValid(*atomStoragePtr_,
                                     domainPtr_->communicator(), true);
      }
      #endif // ifdef DDMD_EXCHANGER_DEBUG
      #endif // ifdef UTIL_DEBUG

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

            if (gridFlags_[i]) {

               // Pack ghost positions for sending
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::UPDATE);
               size = sendArray_(i, j).size();
               for (k = 0; k < size; ++k) {
                  atomPtr = &sendArray_(i, j)[k];
                  atomPtr->packUpdate(*bufferPtr_);
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
                  atomPtr->unpackUpdate(*bufferPtr_);
                  if (shift) {
                     boundaryPtr_->applyShift(atomPtr->position(), i, shift);
                  }
               }
               bufferPtr_->endRecvBlock();
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

            if (gridFlags_[i]) {

               // Pack ghost forces for sending
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::FORCE);
               size = recvArray_(i, j).size();
               for (k = 0; k < size; ++k) {
                  atomPtr = &recvArray_(i, j)[k];
                  atomPtr->packForce(*bufferPtr_);
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
                  atomPtr->unpackForce(*bufferPtr_);
               }
               bufferPtr_->endRecvBlock();
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
