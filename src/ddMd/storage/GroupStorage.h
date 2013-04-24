#ifndef DDMD_GROUP_STORAGE_H
#define DDMD_GROUP_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
#include <util/format/Int.h>
#include <util/mpi/MpiLoader.h>  
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
      * Set parameters, allocate memory and initialize.
      *
      * Call on all processors. Call either this or initialize, but not 
      * both. The initialize function is provided for unit testing.
      *
      * \param capacity      max number of groups owned by processor.
      * \param totalCapacity max number of groups on all processors.
      */
      void initialize(int capacity, int totalCapacity);

      /**
      * Read parameters, allocate memory and initialize.
      *
      * Call on all processors. Call either this or initialize, but not both.
      *
      * Parameter file format:
      *  - capacity      [int]  max number of groups owned by processor.
      *  - totalcapacity [int]  max number of groups on all processors.
      *
      * \param in input stream from which to read parameters.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * Call on all processors.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * Call only on ioProcessor (master).
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
      /// \name Group Management
      //@{

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
      * Remove all groups.
      */
      void clearGroups(); 

      //@}
      /// \name Iterator Interface
      //@{
      
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

      //@}
      /// \name Accessors
      //@{

      /**
      * Find local Group<N> indexed by global id.
      * 
      * \return pointer to Group<N> object, or null pointer if absent.
      */
      Group<N>* find(int id) const;

      /**
      * Return current number of groups on this processor.
      */
      int size() const;

      /**
      * Return capacity for groups on this processor.
      */
      int capacity() const;

      /**
      * Return maximum allowable number of groups on all processors.
      *
      * Note: Group ids must be in range 0, ..., totalCapacity-1
      */
      int totalCapacity() const;

      /**
      * Compute and store the number of distinct groups on all processors.
      *
      * This is an MPI reduce operation. The correct result is stored only
      * on the rank 0 processor. 
      *
      * Algorithm: For purposes of counting, each group is assigned to the
      * processor that owns its first atom (index 0), and then values from
      * different processors are summed and stored on the master.
      *
      * \param  communicator MPI communicator for this system.
      * \return on master node, return total number of groups.
      */
      #ifdef UTIL_MPI
      void computeNTotal(MPI::Intracomm& communicator);
      #else
      void computeNTotal();
      #endif
   
      /**
      * Return total number of distinct groups on all processors.
      *
      * This function should be called only on the master processors, after
      * a subsequent call to computeNTotal().
      */
      int nTotal() const;

      /**
      *  Mark nTotal as unknown.
      */
      void unsetNTotal();

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
      * \param atomStorage  associated AtomStorage object
      * \param hasGhosts    true if the atomStorage has ghosts, false otherwise
      * \param communicator domain communicator 
      */
      #ifdef UTIL_MPI
      bool isValid(AtomStorage& atomStorage, MPI::Intracomm& communicator, 
                   bool hasGhosts);
      #else
      bool isValid(AtomStorage& atomStorage, bool hasGhosts);
      #endif

      //@}
      /// \name Statistics
      //@{
 
      /**
      * Compute statistics (reduce from all processors).
      * 
      * Call on all processors.
      */
      #ifdef UTIL_MPI
      virtual void computeStatistics(MPI::Intracomm& communicator);
      #else
      virtual void computeStatistics();
      #endif

      /**
      * Clear statistical accumulators (call on all processors).
      */
      void clearStatistics();

      /**
      * Output statistics.
      *
      * Call on master, after calling computeStatistics on all procs.
      *
      * \param out   output stream
      */
      void outputStatistics(std::ostream& out);

      /**
      * Get the maximum number of primary atoms encountered thus far.
      *
      * Call only on master.
      */
      int maxNGroup() const;

      //@}

   private:

      // Array that holds all available group objects.
      DArray< Group<N> >     groups_;

      // Set of pointers to local groups.
      ArraySet< Group<N> >   groupSet_;

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

      /// Maximum of nAtom1_ on this proc since stats cleared.
      int  maxNGroupLocal_;     
   
      /// Maximum of nAtom1_ on all procs (defined only on master).
      Setable<int>  maxNGroup_;     
      
      // Total number of distinct groups on all processors.
      Setable<int> nTotal_;

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
   {  return nTotal_.value(); }

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
      nTotal_(0)
   {}
 
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

}
#endif
