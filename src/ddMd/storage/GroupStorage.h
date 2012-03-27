#ifndef DDMD_GROUP_STORAGE_H
#define DDMD_GROUP_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <ddMd/chemistry/Group.h>        // member template parameter
#include <util/containers/DArray.h>      // member template
#include <util/containers/ArraySet.h>    // member template
#include <util/containers/ArrayStack.h>  // member template

#include "AtomStorage.h"
#include "GroupIterator.h"
#include "ConstGroupIterator.h"
#include <util/global.h>


namespace DdMd
{

   using namespace Util;

   /**
   * A container for all the Group<N> objects on this processor.
   *
   * \ingroup DdMd_Storage_Module
   */
   template <int N>
   class GroupStorage : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      GroupStorage();

      /**
      * Destructor.
      */
      ~GroupStorage();

      /**
      * Read parameters, allocate memory and initialize.
      *
      * Parameters (file format):
      *  - capacity      [int]  max number of groups owned by processor.
      *  - totalcapacity [int]  max number of groups on all processors.
      */
      virtual void readParam(std::istream& in);

      /**
      * Set parameters, allocate memory and initialize.
      *
      * \param capacity      max number of groups owned by processor.
      * \param totalcapacity max number of groups on all processors.
      */
      void setParam(int capacity, int totalCapacity);

      /**
      * Set value for total number of distinct groups on all processors.
      *
      * This value is used by the isValid() method, and should be set only 
      * on the master processor.
      */
      void setNTotal(int nTotal);

      // Mutators

      /**
      * Returns pointer to an address available for a new Group.
      *
      * This method returns the address of an empty Group<N> object that
      * can be used for a new local Group<N>. The pointer that it returns
      * is popped off the reservoir of unused objects. This method does 
      * not modify the object id or add the Group to the set of local 
      * Group<N> objects, and so must be followed by either a call to 
      * add(), which adds the object to the set, or returnPtr(), which
      * pushes the new pointer back onto the reservoir.
      *
      * \return address available for a new group.
      */
      Group<N>* newPtr();

      /**
      * Reverses the action of newPtr.
      *
      * This method pushes the pointer returned by a previous call to
      * newPtr back onto the reservoir, without adding it to the local
      * group set.
      */
      void returnPtr();

      /**
      * Adds the most recent new Group to the set of local groups.
      *
      * This method adds the Group that is pointed at by the pointer
      * returned by the most recent call to newPtr(). This method
      * thus completes the process initiated by newPtr().
      *
      * Usage:
      * \code
      * GroupStorage<N> storage;
      *
      * Group<N>* ptr = storage.newPtr();
      * ptr->setId(id);
      * storage.add();
      *
      * \endcode
      * The Group id must be set before calling add(), but other
      * properties may be set either before or after calling add().
      *
      * Preconditions:
      * 1) newPtr() must have been called since the last add().
      * 2) A valid global group id for this pointer must have been set.
      * 3) A Group with the specified id must not already be in the set.
      *
      * \return pointer to newly added Group.
      */
      Group<N>* add();

      /**
      * Add a new Group with a specified id.
      *
      * Adds a new Group<N> to the set owned by this method, with a
      * specified global id, and returns a pointer to the new Group.
      * This method calls newPtr(), and sets the id of the new group,
      * and calls add().  Other member variables must be set after
      * calling this method, using the returned pointer.
      * 
      * \param id  global id for the new Group.
      * \return pointer to the new Group.
      */
      Group<N>* add(int id); 

      /**
      * Remove a specific Group.
      *
      * \throw Exception if groupPtr is not address of a Group.
      *
      * \param groupPtr pointer to the group to be removed
      */
      void remove(Group<N>* groupPtr); 

      /**
      * Add to set of incomplete groups.
      *
      * \throw Exception if groupPtr is not address of a local Group.
      *
      * \param groupPtr pointer to the group to be added to set.
      */
      void addIncomplete(Group<N>* groupPtr); 

      /**
      * Remove from set of incomplete groups.
      *
      * \throw Exception if groupPtr is not address of a Group.
      *
      * \param groupPtr pointer to the group to be removed
      */
      void removeIncomplete(Group<N>* groupPtr); 

      /**
      * Clear the set of incomplete groups.
      */
      void clearIncompleteSet(); 

      /**
      * Compute and store the number of distinct groups on all processors.
      *
      * This is an MPI reduce operation. The correct result is stored and
      * returned only on the rank 0 processor. On other processors, the
      * method stores a null value of -1.
      *
      * \param  communicator MPI communicator for this system.
      * \return on master node, return total number of groups.
      */
      #ifdef UTIL_MPI
      void computeNTotal(MPI::Intracomm& communicator);
      #else
      void computeNTotal();
      #endif
   
      // Accessors

      /**
      * Find local Group<N> indexed by global id.
      * 
      * \return pointer to Group<N> object, or null pointer if absent.
      */
      Group<N>* find(int id) const;

      /**
      * Set iterator to beginning of the set of groups.
      *
      * \param iterator iterator for all groups.
      */
      void begin(GroupIterator<N>& iterator);
 
      /**
      * Set iterator to beginning of the set of groups.
      *
      * \param iterator iterator for all groups.
      */
      void begin(ConstGroupIterator<N>& iterator) const;
 
      /**
      * Set iterator to beginning of the set of incomplete groups.
      *
      * \param iterator iterator for all incomplete groups.
      */
      void begin(IncompleteGroupIterator<N>& iterator);
 
      /**
      * Set const iterator to beginning of the set of incomplete groups.
      *
      * \param iterator iterator for all incomplete groups.
      */
      void begin(ConstIncompleteGroupIterator<N>& iterator) const;
 
      /**
      * Return current number of groups on this processor.
      */
      int size() const;

      /**
      * Return capacity for groups on this processor.
      */
      int capacity() const;

      /**
      * Return maximum number of groups on all processors.
      *
      * Group ids are labelled from 0, ..., totalCapacity-1
      */
      int totalCapacity() const;

      /**
      * Return total number of distinct groups on all processors.
      */
      int nTotal() const;

      /**
      * Return true if the container is valid, or throw an Exception.
      */
      bool isValid();

      /**
      * Return true if the container is valid, or throw an Exception.
      *
      * Calls overloaded isValid() method, then checks consistency of atom 
      * pointers with those in an asociated AtomStorage. If hasGhosts is
      * false, the method requires that no group contain a pointer to a ghost 
      * atom. If hasGhosts is true, requires that every Group be complete.
      *
      * \param atomStorage associated AtomStorage object.
      * \param hasGhosts   true if the atomStorage has ghosts, false otherwise.
      */
      #ifdef UTIL_MPI
      bool isValid(AtomStorage& atomStorage, MPI::Intracomm& communicator, 
                   bool hasGhosts);
      #else
      bool isValid(AtomStorage& atomStorage, bool hasGhosts);
      #endif

   private:

      // Array that holds all available group objects.
      DArray< Group<N> >     groups_;

      // Set of pointers to local groups.
      ArraySet< Group<N> >   groupSet_;

      // Set of pointers to incomplete local groups.
      ArraySet< Group<N> >   incompleteGroupSet_;

      // Stack of pointers to unused local Group objects.
      ArrayStack< Group<N> > reservoir_;

      // Array of pointers to groups, indexed by Id.
      // Elements corresponding to absent groups hold null pointers.
      DArray< Group<N>* >    groupPtrs_;

      // Pointer to space for a new local Group
      Group<N>* newPtr_;

      // Capacity for local groups on this processor.
      int capacity_;

      // Maximum number of groups on all processors, maximum id + 1
      int totalCapacity_;

      // Total number of distinct groups on all processors.
      int nTotal_;

      /*
      * Allocate and initialize all private containers.
      */
      void allocate();
    
   };

   // Inline method definitions

   template <int N>
   inline int GroupStorage<N>::size() const
   {  return groupSet_.size(); }

   template <int N>
   inline int GroupStorage<N>::capacity() const
   {  return capacity_; }

   template <int N>
   inline int GroupStorage<N>::totalCapacity() const
   {  return totalCapacity_; }

   template <int N>
   inline int GroupStorage<N>::nTotal() const
   {  return nTotal_; }

   // Non-inline method templates.

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
      nTotal_(-1)
   {}
 
   /*
   * Destructor.
   */
   template <int N>
   GroupStorage<N>::~GroupStorage()
   {}

   /*
   * Read parameters and allocate memory.
   */
   template <int N>
   void GroupStorage<N>::readParam(std::istream& in)
   {
      read<int>(in, "capacity", capacity_);
      read<int>(in, "totalCapacity", totalCapacity_);
      allocate();
   }

   /*
   * Set parameters and allocate memory.
   */
   template <int N>
   void GroupStorage<N>::setParam(int capacity, int totalCapacity)
   {
      capacity_      = capacity;
      totalCapacity_ = totalCapacity;
      allocate();
   }

   /*
   * Set total number of distinct groups on all processors.
   */
   template <int N>
   void GroupStorage<N>::setNTotal(int nTotal)
   {  nTotal_ = nTotal; }

   /*
   * Allocate and initialize all containers (private).
   */
   template <int N>
   void GroupStorage<N>::allocate()
   {
      groups_.allocate(capacity_);
      reservoir_.allocate(capacity_);
      groupSet_.allocate(groups_);
      incompleteGroupSet_.allocate(groups_);
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
   * Add to set of incomplete groups.
   */
   template <int N>
   void GroupStorage<N>::addIncomplete(Group<N>* groupPtr)
   {
      int id = groupPtr->id();
      if (groupPtrs_[id] == 0) {
         UTIL_THROW("Ptr does not point to local group");
      }
      incompleteGroupSet_.append(*groupPtr);
   }

   /*
   * Remove from set of incomplete groups.
   */
   template <int N>
   void GroupStorage<N>::removeIncomplete(Group<N>* groupPtr)
   {
      int id = groupPtr->id();
      if (groupPtrs_[id] == 0) {
         UTIL_THROW("Ptr does not point to local group");
      }
      incompleteGroupSet_.remove(*groupPtr);
   }

   /*
   * Clear incomplete group set (remove all).
   */
   template <int N>
   void GroupStorage<N>::clearIncompleteSet()
   {  incompleteGroupSet_.clear(); }

   // Accessors

   /*
   * Return pointer to Group with specified id.
   */
   template <int N>
   Group<N>* GroupStorage<N>::find(int id) const
   {  return groupPtrs_[id]; }

   /*
   * Set iterator to beginning of the set of local groups.
   */
   template <int N>
   void GroupStorage<N>::begin(GroupIterator<N>& iterator)
   {  groupSet_.begin(iterator); }
 
   /*
   * Set const iterator to beginning of the set of local groups.
   */
   template <int N>
   void GroupStorage<N>::begin(ConstGroupIterator<N>& iterator) const
   {  groupSet_.begin(iterator); }
 
   /*
   * Set iterator to beginning of the set of incomplete local groups.
   */
   template <int N>
   void GroupStorage<N>::begin(IncompleteGroupIterator<N>& iterator)
   {  incompleteGroupSet_.begin(iterator); }
 
   /*
   * Set const iterator to beginning of the set of incomplete local groups.
   */
   template <int N>
   void GroupStorage<N>::begin(ConstIncompleteGroupIterator<N>& iterator) const
   {  incompleteGroupSet_.begin(iterator); }
 
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
      for (begin(iter); !iter.atEnd(); ++iter) {
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
   * Compute, store and return total number of atoms on all processors.
   */
   template <int N>
   #ifdef UTIL_MPI
   void GroupStorage<N>::computeNTotal(MPI::Intracomm& communicator)
   #else
   void GroupStorage<N>::computeNTotal()
   #endif
   {
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
      #ifdef UTIL_MPI
      communicator.Reduce(&nLocal, &nTotal_, 1, 
                          MPI::INT, MPI::SUM, 0);
      if (communicator.Get_rank() !=0) {
         nTotal_ = -1;
      }
      #else
      nTotal_ = nLocal;
      #endif
   }

   /*
   * Check validity of all groups on this processor.
   */
   template <int N>
   #ifdef UTIL_MPI
   bool GroupStorage<N>::isValid(AtomStorage& atomStorage, MPI::Intracomm& communicator,
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

      #ifdef UTIL_MPI
      // Count & return number of local atoms in groups on all processors.
      int nAtomGroupTotal;
      const int source = 0;
      communicator.Reduce(&nAtomGroup, &nAtomGroupTotal, 1, 
                          MPI::INT, MPI::SUM, source);
      if (communicator.Get_rank() == source) {
         if (nTotal_ < 0) {
            UTIL_THROW("nTotal not set before isValid");
         }
         if (nAtomGroupTotal != nTotal_*N) {
            std::cout << "nAtomGroupTotal = " << nAtomGroupTotal << std::endl;
            std::cout << "nTotal*N        = " << nTotal_*N  << std::endl;
            UTIL_THROW("Discrepancy in number of local atoms in Group objects");
         }
      }
      #endif

      return true;
   }
}
#endif
