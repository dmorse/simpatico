#ifndef DDMD_GROUP_STORAGE_TPP
#define DDMD_GROUP_STORAGE_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "GroupStorage.h"
#include "AtomStorage.h"
#include "AtomMap.h"
#include <util/format/Int.h>
#include <util/mpi/MpiLoader.h>  

//#define DDMD_GROUP_STORAGE_DEBUG

namespace DdMd
{
 
   using namespace Util;

   /*
   * Default constructor.
   */
   template <int N>
   GroupStorage<N>::GroupStorage()
    : groups_(),
      groupSet_(),
      reservoir_(),
      newPtr_(0),
      nGroupDistinct_(0),
      capacity_(0),
      totalCapacity_(0),
      nTotal_(0)
   {  emptyGroups_.reserve(128); }
 
   /*
   * Destructor.
   */
   template <int N>
   GroupStorage<N>::~GroupStorage()
   {}

   /*
   * Set parameters and allocate memory.
   */
   template <int N>
   void GroupStorage<N>::initialize(int capacity, int totalCapacity)
   {
      capacity_  = capacity;
      totalCapacity_ = totalCapacity;
      allocate();
   }

   /*
   * Read parameters and allocate memory.
   */
   template <int N>
   void GroupStorage<N>::readParameters(std::istream& in)
   {
      read<int>(in, "capacity", capacity_);
      read<int>(in, "totalCapacity", totalCapacity_);
      allocate();
   }

   /*
   * Load parameters from input archive and allocate memory.
   */
   template <int N>
   void GroupStorage<N>::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<int>(ar, "capacity", capacity_);
      loadParameter<int>(ar, "totalCapacity", totalCapacity_);
      allocate();

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(maxNGroupLocal_);
      maxNGroup_.set(maxNGroupLocal_);
   }

   /*
   * Save parameters to output archive.
   */
   template <int N>
   void GroupStorage<N>::save(Serializable::OArchive& ar)
   {
      ar & capacity_;
      ar & totalCapacity_;
      int max = maxNGroup_.value();
      ar & max;
   }

   /*
   * Allocate and initialize all containers (private).
   */
   template <int N>
   void GroupStorage<N>::allocate()
   {
      groups_.allocate(capacity_);
      reservoir_.allocate(capacity_);
      groupSet_.allocate(groups_);
      groupPtrs_.allocate(totalCapacity_);

      // Push all groups onto reservoir stack, in reverse order.
      for (int i = capacity_ - 1; i >=0; --i) {
          reservoir_.push(groups_[i]);
      }

      // Nullify all pointers in groupPtrs_ array.
      for (int i = 0; i < totalCapacity_; ++i) {
         groupPtrs_[i] = 0;
      }
   }

   // Local group mutators

   /*
   * Returns address for a new local Group.
   */ 
   template <int N>
   Group<N>* GroupStorage<N>::newPtr()
   {
      // Precondition
      if (newPtr_ != 0) 
         UTIL_THROW("Unregistered newPtr_ still active");
      newPtr_ = &reservoir_.pop();
      newPtr_->clear();
      return newPtr_;
   }

   /*
   * Pushes unused pointer back onto reservoir.
   */ 
   template <int N>
   void GroupStorage<N>::returnPtr()
   {
      // Preconditions
      if (newPtr_ == 0) 
         UTIL_THROW("No active newPtr_");
      newPtr_->setId(-1);
      reservoir_.push(*newPtr_);
      newPtr_ = 0;
   }

   /*
   * Register new local Group in internal data structures.
   */ 
   template <int N>
   Group<N>* GroupStorage<N>::add()
   {

      // Preconditions
      if (newPtr_ == 0) {
         UTIL_THROW("No active newPtr_");
      }
      int groupId = newPtr_->id();
      if (groupId < 0 || groupId >= totalCapacity_) {
         Log::file() << "groupId = " << groupId << std::endl;
         UTIL_THROW("Invalid group id");
      }
      if (groupPtrs_[groupId] != 0) {
         UTIL_THROW("Group with specified id is already present");
      }

      // Add Group<N> object to container
      groupSet_.append(*newPtr_);
      groupPtrs_[groupId] = newPtr_;
      ++nGroupDistinct_;

      // Release newPtr_ for reuse.
      Group<N>* ptr = newPtr_;
      newPtr_ = 0;

      // Check maximum.
      if (groupSet_.size() > maxNGroupLocal_) {
         maxNGroupLocal_ = groupSet_.size();
      }

      return ptr;
   }

   /*
   * Add a new Group with a specified id, return pointer to new Group.
   */
   template <int N>
   Group<N>* GroupStorage<N>::add(int id)
   {
      Group<N>* ptr = newPtr();
      ptr->setId(id);
      add();
      return ptr;
   }

   /*
   * Remove a specific local Group.
   */
   template <int N>
   void GroupStorage<N>::remove(Group<N>* groupPtr)
   {
      int groupId = groupPtr->id();
      if (groupId < 0 || groupId >= totalCapacity_) {
         Log::file() << "Group id = " << groupId << std::endl;
         UTIL_THROW("Invalid group id, out of range");
      } else if (groupPtrs_[groupId] == 0) {
         UTIL_THROW("Group does not exist on this processor");
      }
      reservoir_.push(*groupPtr);
      groupSet_.remove(*groupPtr);
      groupPtrs_[groupId] = 0;
      --nGroupDistinct_;
      groupPtr->setId(-1);
   }

   /*
   * Remove all groups.
   */
   template <int N>
   void GroupStorage<N>::clearGroups()
   {
      Group<N>* groupPtr;
      int  groupId;
      while (groupSet_.size() > 0) {
         groupPtr = &groupSet_.pop();
         groupId = groupPtr->id();
         groupPtrs_[groupId] = 0;
         groupPtr->setId(-1);
         reservoir_.push(*groupPtr);
      }
      if (groupSet_.size() != 0) {
         UTIL_THROW("Nonzero groupSet size at end of clearGhosts");
      }
      nGroupDistinct_ = 0;
      ghosts_.clear();
   }

   // Accessors

   /**
   * Compute and store total number of atoms on all processors.
   */
   template <int N>
   #ifdef UTIL_MPI
   void GroupStorage<N>::computeNTotal(MPI::Intracomm& communicator)
   #else
   void GroupStorage<N>::computeNTotal()
   #endif
   {
      // If nTotal is already known, return and do nothing.
      if (nTotal_.isSet()) return;

      // Loop over all groups on this processor (including ghosts).
      // Increment nLocal only if atom 0 is local on this processor
      GroupIterator<N> iterator;
      Atom* atomPtr;
      int nLocal = 0;
      begin(iterator);
      for ( ; iterator.notEnd(); ++iterator) {
         atomPtr = iterator->atomPtr(0);
         if (atomPtr) {
            if (!atomPtr->isGhost()) {
               ++nLocal;
            }
         }
      }

      // Reduce data on all processors and set nTotal_ on master.
      int sum;
      #ifdef UTIL_MPI
      communicator.Reduce(&nLocal, &sum, 1, MPI::INT, MPI::SUM, 0);
      if (communicator.Get_rank() !=0) {
         sum = -1;
      }
      nTotal_.set(sum);
      #else
      nTotal_.set(nLocal);
      #endif
   }

   /*
   * Compute memory usage statistics (call on all processors).
   */
   template <int N>
   #ifdef UTIL_MPI
   void GroupStorage<N>::computeStatistics(MPI::Intracomm& communicator)
   #else
   void GroupStorage<N>::computeStatistics()
   #endif
   { 
      #ifdef UTIL_MPI
      int maxNGroupGlobal;
      communicator.Allreduce(&maxNGroupLocal_, &maxNGroupGlobal, 1, 
                             MPI::INT, MPI::MAX);
      maxNGroup_.set(maxNGroupGlobal);
      maxNGroupLocal_ = maxNGroupGlobal;
      #else
      maxNGroup_.set(maxNGroupLocal_);
      #endif
   }

   /*
   * Clear all statistics.
   */
   template <int N>
   void GroupStorage<N>::clearStatistics() 
   {
      maxNGroupLocal_ = 0;
      maxNGroup_.unset();
   }

   /*
   * Output statistics.
   */
   template <int N>
   void GroupStorage<N>::outputStatistics(std::ostream& out)
   {

      out << std::endl;
      out << "GroupStorage<" << N << ">" << std::endl;
      out << "NGroup: max, capacity    " 
                  << Int(maxNGroup_.value(), 10)
                  << Int(capacity_, 10)
                  << std::endl;
   }

   /*
   * Check validity of this GroupStorage (minimal checks).
   *
   * Returns true if all is ok, or throws an Exception.
   */
   template <int N>
   bool GroupStorage<N>::isValid()
   {
      if (size() + reservoir_.size() != capacity_) {
         UTIL_THROW("nGroup + reservoir size != local capacity"); 
      }

      // Check consistency of group ids and indexing in groupPtrs_
      Group<N>* ptr;
      int       i, j;
      j = 0;
      for (i = 0; i < totalCapacity_ ; ++i) {
         ptr = groupPtrs_[i];
         if (ptr != 0) {
            ++j;
            if (ptr->id() != i) {
               UTIL_THROW("ptr->id() != i"); 
            }
         }
      }

      // Count local groups
      GroupIterator<N> iter;
      j = 0;
      for (begin(iter); iter.notEnd(); ++iter) {
         ++j;
         ptr = find(iter->id());
         if (ptr == 0) {
            UTIL_THROW("Unable to find local group returned by iterator"); 
         }
         if (ghosts_.size() == 0) {
            if (ptr != iter.get()) {
               UTIL_THROW("Inconsistent find(iter->id()"); 
            }
         }
      }
      if (j != size()) {
         UTIL_THROW("Number from iterator != size()"); 
      }

      return true;
   }

   /*
   * Check validity of all groups on this processor.
   */
   template <int N>
   bool
   #ifdef UTIL_MPI
   GroupStorage<N>::isValid(const AtomStorage& atomStorage, 
                            MPI::Intracomm& communicator, bool hasGhosts)
   #else
   GroupStorage<N>::isValid(const AtomStorage& atomStorage, bool hasGhosts)
   #endif
   {
      // Call simple isValid() function to check only GroupStorage data.
      isValid();

      int i;
      int atomId;
      int nAtom;  // # of nonnull pointers to atoms in one group
      int nLocal; // # of local atoms in one group
      int nGroupLocal = 0; // # of local atoms in all groups on processor
      Atom* atomPtr;
      Atom* findPtr;
      ConstGroupIterator<N> groupIter;

      // Loop over groups.
      const AtomMap& atomMap = atomStorage.map();
      for (begin(groupIter); groupIter.notEnd(); ++groupIter) {
         nAtom = 0;
         nLocal = 0;
         for (i = 0; i < N; ++i) {
            atomId  = groupIter->atomId(i);
            if (atomId < 0 || atomId >= atomStorage.totalAtomCapacity()) {
               UTIL_THROW("Invalid atom id in Group");
            }
            findPtr = atomMap.find(atomId);
            atomPtr = groupIter->atomPtr(i);
            if (atomPtr) {
               ++nAtom;
               if (atomPtr->id() != atomId) {
                  UTIL_THROW("Inconsistent id for atom pointer in Group");
               }
               if (findPtr == 0) {
                  UTIL_THROW("Pointer to atom that is not in the map");
               }
               if (!atomPtr->isGhost()) {
                  ++nLocal; 
                  if (atomPtr != findPtr) {
                     UTIL_THROW("Inconsistent pointer to local atom");
                  }
               } else { // if atomPtr points to a ghost
                  if (!hasGhosts) {
                     UTIL_THROW("Unexpected pointer to ghost atom");
                  }
               }
            } else { // if atomPtr is null
               if (hasGhosts) {
                  UTIL_THROW("Unexpected null pointer in group");
               }
               if (findPtr) {
                  if (!findPtr->isGhost()) {
                     UTIL_THROW("Null atomPtr but local atom in map");
                  }
               }
            }
         }
         if (nAtom == 0) {
            UTIL_THROW("Empty group");
         }
         if (hasGhosts && nAtom < N) {
            UTIL_THROW("Incomplete group");
         }
         nGroupLocal += nLocal;
      }

      // Count number of distinct groups.
      #ifdef UTIL_MPI
      unsetNTotal();
      computeNTotal(communicator);
      #endif

      #ifdef UTIL_MPI
      // Compute number of local atoms in groups on all processors.
      int sum; // Total number of local atoms
      const int source = 0; // source processor id for reduce
      communicator.Reduce(&nGroupLocal, &sum, 1, 
                          MPI::INT, MPI::SUM, source);
      if (communicator.Get_rank() == source) {
         if (!nTotal_.isSet()) {
            UTIL_THROW("nTotal not set");
         }
         if (sum != N*nTotal()) {
            Log::file() << "# of local atoms = " << sum << std::endl;
            Log::file() << "nTotal*N         = " 
                        << N*nTotal() << std::endl;
            UTIL_THROW("Wrong number of local atoms in Groups");
         }
      }
      #endif

      return true;
   }

   /*
   * Check validity, after atom exchange but before ghost exchange.
   *
   * This function may only be called after exchange of atoms and groups,
   * but before exchange and creation of ghosts.
   */
   template <int N>
   bool
   #ifdef UTIL_MPI
   GroupStorage<N>::isValid(const AtomStorage& atomStorage, 
                            MPI::Intracomm& communicator)
   {  return isValid(atomStorage, communicator, false); }
   #else
   isValid(const AtomStorage& atomStorage)
   {  return isValid(atomStorage, false); }
   #endif

   /*
   * Check validity, after ghost exchange is complete.
   */
   template <int N>
   bool
   #ifdef UTIL_MPI
   GroupStorage<N>::isValid(const AtomStorage& atomStorage, 
                            const Boundary& boundary,
                            MPI::Intracomm& communicator)
   #else
   GroupStorage<N>::isValid(const AtomStorage& atomStorage, 
                            const Boundary& boundary);
   #endif
   {
      assert(!atomStorage.isCartesian());

      bool hasGhosts = true;
      #ifdef UTIL_MPI
      isValid(atomStorage, communicator, hasGhosts);
      #else
      isValid(atomStorage, hasGhosts);
      #endif

      // Check that all groups are spatially compact.
      ConstGroupIterator<N> groupIter;
      for (begin(groupIter); groupIter.notEnd(); ++groupIter) {
         if (!groupIter->isCompactGen(boundary)) {
            UTIL_THROW("Noncompact group");
         }
      }

      return true;
   }

   /*
   *  Mark nTotal as unknown.
   */
   template <int N>
   void GroupStorage<N>::unsetNTotal()
   {  nTotal_.unset(); }

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
   * "down" (j=0) and "up" (j=1) in each direction). This requires information
   * about positions of ghost as well as local atoms. For each boundary of 
   * the domain, identify atoms whose positions are "inside" and "outside".
   * Count ghost atoms very near the boundary as both inside and outside,
   * for safety. If a group has atoms both inside and outside a domain 
   * boundary, it is marked for sending in the associated communication 
   * step. 
   *
   * After calculating a ghost communication plan for each group, clear 
   * the pointers to all ghost atoms in the group. The exchangeAtoms 
   * method will clear the actual ghost atoms from the AtomStorage.
   */
   template <int N> void 
   GroupStorage<N>::beginAtomExchange(FMatrix<double, Dimension, 2>& bound, 
                                      FMatrix<double, Dimension, 2>& inner, 
                                      FMatrix<double, Dimension, 2>& outer, 
                                      IntVector& gridFlags, 
                                      const AtomMap& map)
   {
      double coordinate;
      GroupIterator<N> groupIter;
      Group<N>* groupPtr;
      Atom* atomPtr;
      int nIn, nOut, i, j, k;
      bool isComplete;
      bool choose;

      // Clear ghost groups (if any)
      for (j = 0; j < ghosts_.size(); ++j) {
         groupPtr = &ghosts_[j];
         groupSet_.remove(*groupPtr);
         reservoir_.push(*groupPtr);
         groupPtr->setId(-1);
         #if 0
         for (k = 0; k < N; ++k) {
            groupPtr->setAtomId(k, -1);
            groupPtr->clearAtomPtr(k);
         }
         #endif
      }
      ghosts_.clear();

      // Loop over groups
      begin(groupIter);
      for ( ; groupIter.notEnd(); ++groupIter) {
         groupIter->plan().clearFlags();

         isComplete = (groupIter->nPtr() == N); 
         if (isComplete) { // if group is complete

            for (i = 0; i < Dimension; ++i) {
               for (j = 0; j < 2; ++j) {
 
                  // Determine if Group may span boundary (i, j)
                  choose = false; 
                  nIn = 0;
                  nOut = 0;
                  // Loop over atoms in group
                  for (k = 0; k < N; ++k) {
                     atomPtr = groupIter->atomPtr(k);
                     assert(atomPtr);
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
                  } // end for k (atoms in group)
                  if (nOut > 0 && nIn > 0) {
                     choose = true;
                  }

                  // If group may span boundary (i,j), set ghost flag (i, j)
                  if (choose) {
                     groupIter->plan().setGhost(i, j);
                  }

               } // end for j = 0, 1
  
               // A complete group may not span both lower (j=0) and upper (j=1) boundaries
               if (groupIter->plan().ghost(i, 0) && groupIter->plan().ghost(i, 1)) {
                  Log::file() << "Direction " << i << std::endl;
                  Log::file() << "Inner / outer (j=0) = " << inner(i,0) 
                              << "  " << outer(i, 0) << std::endl;
                  Log::file() << "Inner / outer (j=1) = " << inner(i,1) 
                              << "  " << outer(i, 1) << std::endl;
                  for (k = 0; k < N; ++k) {
                     atomPtr = groupIter->atomPtr(k);
                     assert(atomPtr);
                     coordinate = atomPtr->position()[i];
                     Log::file() << k << "  " << coordinate;
                     if (atomPtr->isGhost()) {
                        Log::file() << " ghost  ";
                     } else {
                        Log::file() << " local  "
                                    << atomPtr->plan().exchange(i, 0) << "  "
                                    << atomPtr->plan().exchange(i, 1);
                     }
                     Log::file() << std::endl;
                     Log::file() << std::endl;
                  }
                  UTIL_THROW("Group spans both upper and lower boundaries");
               }

               // For directions with one processor, send in both directions or neither
               if (!gridFlags[i]) {
                  if (groupIter->plan().ghost(i,0)) {
                     groupIter->plan().setGhost(i,1);
                  } else 
                  if (groupIter->plan().ghost(i,1)) {
                     groupIter->plan().setGhost(i,0);
                  }
               }

            } // end for i (Cartesian axes)

         } else { // if group is not complete

            // If not complete, mark ghost flag for all directions
            for (i = 0; i < Dimension; ++i) {
               for (j = 0; j < 2; ++j) {
                 groupIter->plan().setGhost(i, j);
               }
            }

         } // if-else (isComplete)

         // Clear pointers to all ghost atoms in this group
         for (k = 0; k < N; ++k) {
            atomPtr = groupIter->atomPtr(k);
            if (atomPtr) {
               if (atomPtr->isGhost()) {
                  groupIter->clearAtomPtr(k);
                  atomPtr = map.find(groupIter->atomId(k));
                  if (atomPtr) {
                     if (!atomPtr->isGhost()) {
                        groupIter->setAtomPtr(k, atomPtr);
                     }
                  }
               }
               #ifdef UTIL_DEBUG
               else {
                  assert(atomPtr == map.find(groupIter->atomId(k)));
               }
               #endif
            }
            #ifdef UTIL_DEBUG
            else {
               assert(0 == map.find(groupIter->atomId(k)));
            }
            #endif
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
      assert(ghosts_.size() == 0);

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
      const AtomMap& atomMap = atomStorage.map();
      int groupId;

      buffer.beginRecvBlock();
      while (buffer.recvSize() > 0) {
         newGroupPtr = newPtr();
         newGroupPtr->unpack(buffer);
         groupId = newGroupPtr->id();
         oldGroupPtr = find(groupId);
         if (oldGroupPtr) {
            returnPtr();
            atomMap.findGroupLocalAtoms(*oldGroupPtr);
         } else {
            add();
            atomMap.findGroupLocalAtoms(*newGroupPtr);
         }

         /*
         * Notes:
         * 1) On receipt a group is added to the storage only if it is not
         * already present. 
         * 2) Local atoms for every received group must be re-identified by
         * calling AtomMap::findGroupLocalAtoms() after each communication
         * in order to identify atoms that have just arrived. Pointers to
         * departing atoms are cleared from groups on the sending processor
         * in the GroupStorage::pack() function. Group pointers to local atoms 
         * must be correct after each buffer exchange because they are used 
         * on subsequent steps to determine which groups must be exchanged 
         * when an atom is exchanged.
         */

      }
      buffer.endRecvBlock();
      assert(buffer.recvSize() == 0);
   }
   #endif // endif ifdef UTIL_MPI

   /*
   * Set ghost communication flags for all atoms in incomplete groups.
   *
   * Precondition: This is called by exchangeAtoms after exchanging atoms 
   * groups between neighboring processors. At this point, there are no 
   * no ghosts atoms.
   *
   * Algorithm: Loop over all Group<N> objects in the group storage. For
   * each group, check if the group is incomplete or not compact. A group
   * is incomplete if it has one or more null pointers, implying that 
   * those atoms are owned by another processor. A group is not compact 
   * if it marked for sending in a direction along which the grid has 
   * only one processor, implying that is spans a periodic boundary.
   * If the group is incomplete or not compact, loop over 6 transfer 
   * directions. For each direction (i,j), if the group is marked for 
   * sending in that direction, set the ghost ghost communication flag 
   * for transfer in that direction for every local atom in the group, 
   * and also add each such atom to sendArray(i, j) for that direction.
   */
   template <int N> void 
   GroupStorage<N>::beginGhostExchange(AtomStorage& atomStorage, 
                          FMatrix<GPArray<Atom>, Dimension, 2>& sendArray,
                          IntVector& gridFlags)
   {
      GroupIterator<N> groupIter;
      Atom* atomPtr;
      Plan* planPtr;
      int i, j, k, nAtom;
      bool choose, isPeriodic;

      // Set isPeriodic=true if there is only one processor in any direction
      isPeriodic = false;
      for (i = 0; i < Dimension; ++i) {
         if (!gridFlags[i]) {
            isPeriodic = true;
         }
      }

      // Loop over groups
      begin(groupIter);
      for ( ; groupIter.notEnd(); ++groupIter) {

         #ifdef UTIL_DEBUG
         #ifdef DDMD_GROUP_STORAGE_DEBUG
         // Validate group
         const AtomMap& atomMap = atomStorage.map();
         int atomId;
         nAtom  = 0;
         for (k = 0; k < N; ++k) {
            atomPtr = groupIter->atomPtr(k);
            atomId  = groupIter->atomId(k);
            if (atomPtr != 0) {
               if (atomPtr->isGhost()) {
                  UTIL_THROW("Pointer to ghost atom in group");
               }
               if (atomPtr != atomMap.find(atomId)) {
                  UTIL_THROW("Inconsistent pointer to local atom in group");
               }
               ++nAtom;
            } else { // if atomPtr == 0
               atomPtr = atomMap.find(atomId);
               if (atomPtr) {
                  if (!atomPtr->isGhost()) {
                     UTIL_THROW("Null pointer in group to existing local atom");
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

         // Choose this group if it is incomplete or not compact 
         choose = false;
         nAtom = groupIter->nPtr();
         if (nAtom < N) { // if incomplete
            choose = true; 
         } else 
         if (isPeriodic) { // If the grid has one processor in any direction
            for (i = 0; i < Dimension; ++i) {
               if (!gridFlags[i]) {
                  if (groupIter->plan().ghost(i, 0)) {
                     choose = true;
                  }
               }
            }
         }

         /*
         * If this group is chosen, for each direction (i,j) in which 
         * group ghost is set, loop over atoms in group. For each atom,
         * if ghost flag (i, j) is not already set, set ghost flag (i,j) 
         * and add the atom to the sendArray.
         */ 
         if (choose) {
            for (i = 0; i < Dimension; ++i) {
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



   /*
   * Find ghost members of bonds after exchanging all ghosts.
   *
   * Generic template, used for N=3 and N=4.
   */
   template <int N>
   void GroupStorage<N>::finishGhostExchange(const AtomMap& map, 
                                             const Boundary& boundary)
   {
      GroupIterator<N> groupIter;
      Atom* aPtr;  // Pointer to atom a (must be local)
      Atom* bPtr;  // Pointer to atom b (may be local or ghost)
      int i, j, k; // ids for atom within group

      // Loop over all groups
      for (begin(groupIter); groupIter.notEnd(); ++groupIter) {

         // Identify root i = index of first local atom
         for (i = 0; i < N; ++i) {
            aPtr = groupIter->atomPtr(i);
            if (aPtr) {
              if (!aPtr->isGhost()) {
                break;
              }
            }
         }
         if (aPtr == 0 || aPtr->isGhost()) {
            UTIL_THROW("No local atom found");
         }
 
         // Iterate up from root 
         if (i < N-1) {
            for (j = i; j < N-1; ++j) {
               k = j + 1;
               map.findNearestImage(groupIter->atomId(k), 
                                    aPtr->position(), boundary, bPtr);
               groupIter->setAtomPtr(k, bPtr);
               aPtr = bPtr;
            } 
         }

         // Iterate down from root 
         if (i > 0) {
            aPtr = groupIter->atomPtr(i);
            for (j = i; j > 0; --j) {
               k = j - 1;
               map.findNearestImage(groupIter->atomId(k), 
                                    aPtr->position(), boundary, bPtr);
               groupIter->setAtomPtr(k, bPtr);
               aPtr = bPtr;
            } 
         }

      } // end loop over groups

   }

   /*
   * Make a new image of a group.
   */
   template <int N>
   void GroupStorage<N>::makeGroupImage(Group<N>& group, 
                                        int rootId, Atom* rootPtr,
                                        const AtomMap& map, 
                                        const Boundary& boundary)
   {

      // Get a new group object and add to group and ghost sets
      Group<N>* newPtr = &reservoir_.pop();
      groupSet_.append(*newPtr);
      ghosts_.append(*newPtr);
      if (groupSet_.size() > maxNGroupLocal_) {
         maxNGroupLocal_ = groupSet_.size();
      }

      // Copy group id, atomIds, & root pointer to new group
      newPtr->setId(group.id());
      for (int k=0; k < N; ++k) {
          newPtr->setAtomId(k, group.atomId(k));
          newPtr->clearAtomPtr(k);
      }
      newPtr->setAtomPtr(rootId, rootPtr);

      Atom* aPtr; // first atom of each bond (already known)
      Atom* bPtr; // second atom of each bond (searched for)
      Atom* cPtr; // unused local image of second atom, if any
      int j, k;   // ids within group of first and second atoms

      // Set pointers, iterating up from rootId 
      if (rootId < N-1) {
         aPtr = rootPtr;
         for (j = rootId; j < N - 1; ++j) {
            k = j + 1;
            cPtr = map.findNearestImage(newPtr->atomId(k), 
                                        aPtr->position(), boundary, bPtr);
            newPtr->setAtomPtr(k, bPtr);
            aPtr = bPtr;
            if (cPtr) {
               makeGroupImage(*newPtr, k, cPtr, map, boundary);
            }
         }
      }

      // Set pointers, iterating down from rootId
      if (rootId > 0) {
         aPtr = rootPtr;
         for (j = rootId; j > 0; --j) {
            k = j - 1;
            cPtr = map.findNearestImage(newPtr->atomId(k), 
                                        aPtr->position(), boundary, bPtr);
            newPtr->setAtomPtr(k, bPtr);
            aPtr = bPtr;
            if (cPtr) {
               makeGroupImage(*newPtr, k, cPtr, map, boundary);
            }
         }
      }

   }

} // namespace DdMd
#endif
