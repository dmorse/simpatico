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
#include <util/format/Int.h>
#include <util/mpi/MpiLoader.h>  

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
         std::cout << "groupId = " << groupId << std::endl;
         UTIL_THROW("Invalid group id");
      }
      if (groupPtrs_[groupId] != 0) {
         UTIL_THROW("Group with specified id is already present");
      }

      // Add Group<N> object to container
      groupSet_.append(*newPtr_);
      groupPtrs_[groupId] = newPtr_;

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
         std::cout << "Group id = " << groupId << std::endl;
         UTIL_THROW("Invalid group id, out of range");
      } else if (groupPtrs_[groupId] == 0) {
         UTIL_THROW("Group does not exist on this processor");
      }
      reservoir_.push(*groupPtr);
      groupSet_.remove(*groupPtr);
      groupPtrs_[groupId] = 0;
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
         UTIL_THROW("Nonzero ghostSet size at end of clearGhosts");
      }
   }

   // Accessors

   /*
   * Check validity of this GroupStorage.
   *
   * Returns true if all is ok, or throws an Exception.
   */
   template <int N>
   bool GroupStorage<N>::isValid()
   {
      
      if (size() + reservoir_.size() != capacity_) 
         UTIL_THROW("nGroup + reservoir size != local capacity"); 

      // Check consitency of pointers to atoms and atom ids
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
         if (ptr == 0)
            UTIL_THROW("Unable to find local group returned by iterator"); 
         if (ptr != iter.get())
            UTIL_THROW("Inconsistent find(iter->id()"); 
      }
      if (j != size())
         UTIL_THROW("Number from iterator != size()"); 

      return true;
   }

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

      // Loop over groups on this processor. 
      // Increment nLocal only if atom 0 is owned by this processor
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
      int nTot;
      #ifdef UTIL_MPI
      communicator.Reduce(&nLocal, &nTot, 1, 
                          MPI::INT, MPI::SUM, 0);
      if (communicator.Get_rank() !=0) {
         nTot = -1;
      }
      nTotal_.set(nTot);
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
   * Check validity of all groups on this processor.
   */
   template <int N>
   #ifdef UTIL_MPI
   bool GroupStorage<N>::isValid(AtomStorage& atomStorage, 
                                 MPI::Intracomm& communicator,
                                 bool hasGhosts)
   #else
   bool GroupStorage<N>::isValid(AtomStorage& atomStorage, bool hasGhosts)
   #endif
   {
      int i;
      int atomId;
      int nAtom;  // Number of local atoms in particular group.
      int nGhost; // Number of local atoms in particular group.
      int nAtomGroup = 0; // Number of local atoms in all groups on processor
      Atom* atomPtr;
      ConstGroupIterator<N> groupIter;

      // Call simpler function that only checks storage data structures.
      isValid();

      // Loop over groups.
      begin(groupIter);
      for ( ; groupIter.notEnd(); ++groupIter) {
         nAtom = 0;
         nGhost = 0;
         for (i = 0; i < N; ++i) {
            atomId  = groupIter->atomId(i);
            if (atomId < 0 || atomId >= atomStorage.totalAtomCapacity()) {
               UTIL_THROW("Invalid atom id in Group");
            }
            atomPtr = groupIter->atomPtr(i);
            if (atomPtr) {
               if (atomPtr != atomStorage.find(atomId)) {
                  UTIL_THROW("Inconsistent non-null atom pointer in Group");
               }
               if (atomPtr->isGhost()) {
                  ++nGhost;
               } else {
                  ++nAtom;
               }
            } else {
               atomPtr = atomStorage.find(atomId);
               if (atomPtr != 0) {
                  if (atomPtr->isGhost()) {
                     if (hasGhosts) {
                          UTIL_THROW("Missing ghost atom");
                     }
                  } else {
                     UTIL_THROW("Missing local atom");
                  }
               }
            }
         }
         if (nAtom == 0) {
            UTIL_THROW("Empty group");
         }
         if (hasGhosts && (nAtom + nGhost) < N) {
            UTIL_THROW("Incomplete group");
         }
         nAtomGroup += nAtom;
      }

      // Count number distinct groups.
      #ifdef UTIL_MPI
      unsetNTotal();
      computeNTotal(communicator);
      #endif

      #ifdef UTIL_MPI
      // Count & return number of local atoms in groups on all processors.
      int nAtomGroupTotal;
      const int source = 0;
      communicator.Reduce(&nAtomGroup, &nAtomGroupTotal, 1, 
                          MPI::INT, MPI::SUM, source);
      if (communicator.Get_rank() == source) {
         if (!nTotal_.isSet()) {
            UTIL_THROW("nTotal not set");
         }
         if (nAtomGroupTotal != N*nTotal()) {
            std::cout << "nAtomGroupTotal = " << nAtomGroupTotal << std::endl;
            std::cout << "nTotal*N        = " << N*nTotal() << std::endl;
            UTIL_THROW("Discrepancy in number of local atoms in Group objects");
         }
      }
      #endif

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
                               FMatrix<GPArray<Atom>, Dimension, 2>& sendArray,
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
                  UTIL_THROW("Error in atom pointer in group");
               }
               if (atomPtr->isGhost()) {
                  UTIL_THROW("Pointer to ghost atom in group");
               } else {
                  ++nAtom;
               }
            } else { // if atomPtr == 0
               atomPtr = atomStorage.find(atomId);
               if (atomPtr) {
                  if (!atomPtr->isGhost()) {
                     UTIL_THROW("Missing pointer to local atom in group");
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
