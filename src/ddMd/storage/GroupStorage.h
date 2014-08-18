#ifndef DDMD_GROUP_STORAGE_H
#define DDMD_GROUP_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/param/ParamComposite.h>   // base class
#include <ddMd/storage/GroupExchanger.h> // base class
#include <ddMd/communicate/GroupDistributor.h>   // member
#include <ddMd/communicate/GroupCollector.h>     // member
#include <ddMd/chemistry/Atom.h>         // member template parameter
#include <ddMd/chemistry/Group.h>        // member template parameter
#include <util/space/IntVector.h>        // member template parameter
#include <util/containers/DArray.h>      // member template
#include <util/containers/GPArray.h>     // member template
#include <util/containers/FMatrix.h>     // member template
#include <util/containers/ArraySet.h>    // member template
#include <util/containers/ArrayStack.h>  // member template

#include "GroupIterator.h"               // inline functions
#include "ConstGroupIterator.h"          // inline functions

namespace DdMd
{

   class Domain;
   class Buffer;
   class AtomStorage;
   using namespace Util;

   /**
   * A container for all the Group<N> objects on this processor.
   *
   * \ingroup DdMd_Storage_Module
   */
   template <int N>
   class GroupStorage : public ParamComposite, public GroupExchanger
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

      /// \name Initialization
      //@{
      
      /**
      * Create associations for distributor and collector members.
      *
      * \param domain Domain object (defines processor grid)
      * \param atomStorage AtomStorage object 
      * \param buffer Buffer object (holds memory for communication)
      */
      void associate(Domain& domain, AtomStorage& atomStorage, 
                     Buffer& buffer);

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
  
      //@}
      /// \name Group Management
      //@{

      /**
      * Returns an address available for addition of a new Group.
      *
      * This function begins a transaction that adds a new Group<N> 
      * object by returning the address of a currently unused Group<N>
      * object. The address of this new object is popped off the 
      * reservoir of unused objects and assigned to an internal "new"
      * pointer (newPtr_) before the function returns. This function 
      * does not modify the object id or add the Group to the set of 
      * local Group<N> objects. It must thus be followed by either a
      * call to add(), which completes the transaction by adding the
      * object to the local set, or by a call to returnPtr(), which
      * reverts the transaction and pushes the new pointer back onto 
      * the reservoir.
      *
      * \return address available for a new group.
      */
      Group<N>* newPtr();

      /**
      * Complete addition of a new Group.
      *
      * This function completes the transaction begun by a previous
      * call to newPtr() by adding the new Group whose address was
      * returned by newPtr() to the set of local groups and nullifying
      * the new Ptr. The invoking function must assign the new group a
      * valid, previously unused global group id before calling the 
      * add() function.
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
      * The group id must be set before calling add(), but other
      * properties may be set either before or after calling add().
      *
      * Preconditions:
      * 1) newPtr() must have been called since the last add().
      * 2) The new Group must have been a assigned a valid global group id.
      * 3) A Group with the specified id may not already be in the set.
      */
      void add();

      /**
      * Add a new Group with a specified id, in a single step.
      *
      * Adds a new Group<N> to the set owned by this function, with a
      * specified global id, and returns a pointer to the new Group.
      * Internally, this function calls newPtr(), and sets the id of 
      * the new group, and calls add().  Other member variables must 
      * be set after calling this function, using the returned pointer.
      * 
      * \param id  global id for the new Group.
      * \return pointer to the new Group.
      */
      Group<N>* add(int id); 

      /**
      * Reverts a transaction begun by the newPtr() function.
      *
      * This function pushes the pointer returned by a previous call 
      * to newPtr() back onto the reservoir of unused objects, without 
      * adding it to the local group set, and the nullifies the "new" 
      * pointer.
      */
      void returnPtr();

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
      *
      * This function resets the GroupStorage to an empty state by
      * pushing the addresses of existing groups onto the reservoir 
      * and then clearing the set of groups.
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
      * This function should be called only on the master processor, 
      * after a previous call to computeNTotal() on all processors.
      */
      int nTotal() const;

      /**
      *  Mark nTotal as unknown.
      *
      *  Call on all processors.
      */
      void unsetNTotal();

      /**
      * Return true if the container is valid, or throw an Exception.
      */
      bool isValid();

      //@}
      /// \name GroupExchanger Interface (Interprocessor Communication)
      //@{
      
      /**
      * Find and mark groups that span boundaries.
      *
      * \param bound  boundaries of domain for this processor
      * \param inner  inner slab boundaries (extended domain of neighbors)
      * \param outer  outer slab boundaries (extended domain of this processor)
      * \param gridFlags  element i is 0 iff gridDimension[i] == 1, 1 otherwise
      */
      virtual
      void markSpanningGroups(FMatrix<double, Dimension, 2>& bound, 
                              FMatrix<double, Dimension, 2>& inner, 
                              FMatrix<double, Dimension, 2>& outer, 
                              IntVector& gridFlags);
   
      #ifdef UTIL_MPI
      /**
      * Pack groups for exchange.
      *
      * Usage: This is called after atom exchange plans are set, 
      * but before atoms makred for exchange in direction i, j
      * have been cleared or removed from the atomStorage.
      *
      * \param i  index of Cartesian communication direction
      * \param j  index for sign of direction
      * \param buffer  Buffer object into which groups are packed
      */
      virtual
      void pack(int i, int j, Buffer& buffer);
   
      /**
      * Unpack groups from buffer and locate available atoms.
      *
      * \param buffer  Buffer object from which groups are unpacked
      * \param atomStorage  AtomStorage used to find atoms pointers
      */
      virtual
      void unpack(Buffer& buffer, AtomStorage& atomStorage);
      #endif // endif ifdef UTIL_MPI
   
      /**
      * Set ghost communication flags for all atoms in incomplete groups.
      *
      * Usage: This is called after exchanging all atoms and groups between 
      * processors, but before exchanging ghosts. At this point, atom
      * ownership is finalized, but there are no ghosts.
      *
      * \param atomStorage AtomStorage object
      * \param sendArray   Matrix of arrays of pointers to ghosts to send
      * \param gridFlags   element i is 0 iff gridDimension[i] == 1, 1 otherwise
      */
      virtual
      void markGhosts(AtomStorage& atomStorage, 
                      FMatrix<GPArray<Atom>, Dimension, 2>&  sendArray,
                      IntVector& gridFlags);
   
      /**
      * Find all ghost members of groups.
      *
      * Usage: This called after all ghosts have been exchanged.
      *
      * \param atomStorage AtomStorage object used to find atom pointers
      */
      virtual
      void findGhosts(AtomStorage& atomStorage);
   
      /**
      * Return true if the container is valid, or throw an Exception.
      *
      * This calls the overloaded isValid() function and then checks 
      * consistency of atom pointers in all groups with those in an
      * asociated AtomStorage. If hasGhosts is false, the function 
      * requires that no group contain a pointer to a ghost atom. If 
      * hasGhosts is true, it requires that every Group be complete.
      *
      * \param atomStorage  associated AtomStorage object
      * \param hasGhosts  true if the atomStorage has ghosts, false otherwise
      * \param communicator  domain communicator 
      */
      #ifdef UTIL_MPI
      virtual
      bool isValid(AtomStorage& atomStorage, MPI::Intracomm& communicator, 
                   bool hasGhosts);
      #else
      virtual
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

      /**
      *  Get the GroupDistributor by reference.
      */
      GroupDistributor<N>& distributor();

      /**
      *  Get the GroupCollector by reference.
      */
      GroupCollector<N>& collector();

   private:

      // Memory pool that holds all available group objects.
      DArray< Group<N> >  groups_;

      // Set of pointers to local groups.
      ArraySet< Group<N> >  groupSet_;

      // Stack of pointers to unused local Group objects.
      ArrayStack< Group<N> >  reservoir_;

      // Array of pointers to groups, indexed by global group Id.
      // Elements corresponding to absent groups hold null pointers.
      DArray< Group<N>* >  groupPtrs_;

      // Array identifying empty groups, marked for later removal 
      GPArray< Group<N> > emptyGroups_;

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

      // Distributor and Collector objects
      GroupDistributor<N> distributor_;
      GroupCollector<N>   collector_;

      /*
      * Allocate and initialize all private containers.
      */
      void allocate();
    
   };

   // Inline member function definitions

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

   /*
   * Return pointer to Group with specified id.
   */
   template <int N>
   inline Group<N>* GroupStorage<N>::find(int id) const
   {  return groupPtrs_[id]; }

   /*
   * Set iterator to beginning of the set of local groups.
   */
   template <int N>
   inline void GroupStorage<N>::begin(GroupIterator<N>& iterator)
   {  groupSet_.begin(iterator); }
 
   /*
   * Set const iterator to beginning of the set of local groups.
   */
   template <int N>
   inline 
   void GroupStorage<N>::begin(ConstGroupIterator<N>& iterator) const
   {  groupSet_.begin(iterator); }

   template <int N>
   inline GroupDistributor<N>& GroupStorage<N>::distributor()
   {  return distributor_; }

   template <int N>
   inline GroupCollector<N>& GroupStorage<N>::collector()
   {  return collector_; }

} // namespace DdMd
#endif // ifndef DDMD_GROUP_STORAGE_H
