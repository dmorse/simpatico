#ifndef DDMD_ATOM_STORAGE_H
#define DDMD_ATOM_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>        // base class

#include <ddMd/chemistry/AtomArray.h>         // member
#include <ddMd/storage/AtomMap.h>             // member
#include <ddMd/communicate/AtomDistributor.h> // member
#include <ddMd/communicate/AtomCollector.h>   // member
#include <ddMd/chemistry/Atom.h>              // member template argument
#include <ddMd/chemistry/Group.h>             // in function templates
#include <simp/boundary/Boundary.h>           // typedef
#include <util/containers/DArray.h>           // member template
#include <util/containers/ArraySet.h>         // member template
#include <util/containers/ArrayStack.h>       // member template
#include <util/misc/Setable.h>                // member template
#include <util/global.h>

class AtomStorageTest;

namespace DdMd
{

   class AtomIterator;
   class ConstAtomIterator;
   class GhostIterator;
   class ConstGhostIterator;

   using namespace Util;
   using namespace Simp;

   /**
   * A container for all the atoms and ghost atoms on this processor.
   *
   * \ingroup DdMd_Storage_Atom_Module
   */
   class AtomStorage : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      AtomStorage();

      /**
      * Destructor.
      */
      ~AtomStorage();

      /**
      * Create associations for distributor and collector.
      *
      * \param domain Domain object (defines processor grid)
      * \param boundary Boundary object (defines periodic unit cell)
      * \param buffer Buffer object (holds memory for communication)
      */
      void associate(Domain& domain, Boundary& boundary, Buffer& buffer);

      /**
      * Set parameters, allocate memory and initialize.
      *
      * Call this or (read|load)Parameters to initialize, but not both.
      *
      * \param atomCapacity      max number of atoms owned by processor.
      * \param ghostCapacity     max number of ghosts on this processor.
      * \param totalAtomCapacity max number of atoms on all processors.
      */
      void initialize(int atomCapacity, int ghostCapacity,
                      int totalAtomCapacity);

      /**
      * Read parameters, allocate memory and initialize.
      *
      * Call either this or initialize(), but not both.
      *
      * Parameters (file format):
      *  - atomCapacity      [int]  max number of atoms owned by processor.
      *  - ghostCapacity     [int]  max number of ghosts on this processor.
      *  - totalatomCapacity [int]  max number of atoms on all processors.
      *
      * \param in input parameter stream.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Set all forces to zero.
      *
      * Sets forces on all local atoms to zero. Also zeros forces for 
      * ghost atoms iff parameter zeroGhosts is true.
      *
      * \param zeroGhosts if true, zero forces on ghost atoms.
      */
      void zeroForces(bool zeroGhosts);
  
      /// \name Serialization (Checkpoint \& Restart)
      //@{

      /**
      * Load internal state from an archive.
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
      /// \name Local Atom Management
      //@{

      /**
      * Returns pointer an address available for a new Atom.
      *
      * This function returns the address of an Atom object that can 
      * be used for a new local Atom. The Atom::clear() function is
      * applied to the new atom before it is returned, so that the
      * id, typeId, isGhost flag, mask, and plan have default values.
      * After this function is called, the storage retains the address
      * of the new atom.  This new atom pointer remains ``active"
      * until a matching call to addNewAtom(), as discussed below.
      *
      * This function does not add the new Atom to the atom set, and so 
      * must be followed by a matching call to addNewAtom() to do so.
      * Usage:
      * \code
      * 
      * AtomStorage storage;
      * Atom*       ptr;
      *
      * ptr = storage.newAtomPtr;
      * ptr->setId(id);
      * // Set other properties of the new Atom.
      * storage.addNewAtom();
      *
      * \endcode
      * The matching call to addNewAtom() adds the new atom to the
      * storage and deactivates the internal new atom pointer.
      *
      * \return address for a new atom
      */
      Atom* newAtomPtr(); 

      /**
      * Finalize addition of the most recent new atom.
      *
      * This function adds the atom that was returned by the most 
      * recent call to newAtomPtr to the atom set. Upon return
      * there is no active new atom pointer. The global atom 
      * id must be set before calling this function, by calling 
      * Atom::setId(int), because the algorithm uses the global
      * id returned by Atom::id().
      */
      void addNewAtom(); 

      /**
      * Add atom with specified global id.
      * 
      * This function adds a new atom to the atom set with a specified
      * atom id, and returns a pointer to the address of the new atom. 
      * It is equivalent to the following, in which storage is an 
      * instance of AtomStorage and ptr is an Atom pointer:
      * \code
      * ptr = storage.newAtomPtr;
      * ptr->setId(id);
      * storage.addNewAtom();
      * \endcode
      * The pointer returned by this function can then be used to set
      * other properties of the new atom. 
      *
      * \param id global index for the new Atom.
      * \return address for new atom.
      */
      Atom* addAtom(int id); 

      /**
      * Remove a specific Atom.
      *
      * \throw Exception if atomPtr is not address of a local Atom.
      *
      * \param atomPtr pointer to the atom to be removed
      */
      void removeAtom(Atom* atomPtr); 

      /**
      * Clear all local atoms.
      */
      void clearAtoms(); 

      /**
      * Return number of local atoms on this procesor (excluding ghosts)
      */
      int nAtom() const;

      /**
      * Return capacity for local atoms on this processor (excluding ghosts).
      */
      int atomCapacity() const;

      //@}
      /// \name Ghost Atom Management
      //@{

      /**
      * Returns pointer an address available for a new ghost Atom.
      *
      * This function returns the address of an Atom object that can 
      * be used for a new ghost Atom. It must be followed by a call
      * to addNewGhost(). Usage:
      * \code
      * 
      * AtomStorage storage;
      * Atom*       ptr;
      *
      * ptr = storage.newGhostPtr;
      * ptr->setId(id);
      * // Set other properties of the new Atom.
      * storage.addNewGhost();
      *
      * \endcode
      *
      * \return address for new ghost atom.
      */
      Atom* newGhostPtr(); 

      /**
      * Register the most recent new ghost atom.
      *
      * This function adds the atom that was returned by the most recent 
      * call to newGhostPtr to the ghost atom set. The global atom 
      * id must be set before calling this function, by calling 
      * Atom::setId(int), because an data structure uses the global
      * id returned by Atom::id(). 
      */
      void addNewGhost(); 

      /**
      * Add ghost atom with specified global id.
      * 
      * This function adds a new atom to the ghost atom set and returns
      * a pointer to the address of the new atom. It is equivalent to
      * the following, in which storage is an AtomStorage object and 
      * ptr is an Atom pointer:
      * \code
      * ptr = storage.newGhostPtr;
      * ptr->setId(id);
      * storage.addNewGhost();
      * \endcode
      *
      * \param id global index for the new Atom.
      * \return address for new ghost atom.
      */
      Atom* addGhost(int id); 

      /**
      * Remove a specific ghost Atom.
      *
      * \throw Exception if atomPtr is not address of a ghost Atom.
      *
      * \param atomPtr pointer to the ghost atom to be removed
      */
      void removeGhost(Atom* atomPtr); 

      /**
      * Clear all ghost atoms.
      */
      void clearGhosts(); 

      /**
      * Return current number of ghost atoms on this processor.
      */
      int nGhost() const;

      /**
      * Return capacity for ghost atoms on this processor.
      */
      int ghostCapacity() const;

      //@}
      /// \name Coordinate Systems
      //@{

      /**
      * Transform positions from Cartesian to generalized coordinates.
      *
      * Transforms position coordinates of all local and ghost atoms.
      *
      * \param boundary periodic boundary conditions
      */
      void transformCartToGen(const Boundary& boundary);
 
      /**
      * Transform positions from generalized to Cartesian coordinates.
      *
      * Transforms position coordinates of all local and ghost atoms.
      *
      * \param boundary periodic boundary conditions
      */
      void transformGenToCart(const Boundary& boundary);

      /**
      * Are atom coordinates Cartesian (true) or generalized (false)?
      */
      bool isCartesian() const;

      //@}
      /// \name Displacement Measurement
      //@{

      /**
      * Record current positions of all local atoms and lock storage.
      * 
      * This function stores positions of local atoms and locks the storage, 
      * prohibiting addition or removal of atoms or ghosts until clearSnapshot 
      * is called.
      */
      void makeSnapshot();

      /**
      * Clear previous snapshot.
      *
      * This function removes the lock imposed by a previous call to
      * makeSnapshot(), allowing changes to atom and ghost sets.
      */
      void clearSnapshot();

      /**
      * Return max-squared displacement since the last snapshot.
      *
      * Note: This is a local operation, and returns only the maximum on this
      * processor. 
      *
      * Throws exception if no valid snapshot is available.
      */
      double maxSqDisplacement();

      //@}
      /// \name Iteration
      //@{

      /**
      * Set iterator to beginning of the set of atoms.
      *
      * \param iterator iterator for all atoms.
      */
      void begin(AtomIterator& iterator);
 
      /**
      * Set iterator to beginning of the set of atoms.
      *
      * \param iterator iterator for all atoms.
      */
      void begin(ConstAtomIterator& iterator) const;
 
      /**
      * Set iterator to beginning of the set of ghost atoms.
      *
      * \param iterator iterator for all ghost atoms.
      */
      void begin(GhostIterator& iterator);

      /**
      * Set iterator to beginning of the set of ghost atoms.
      *
      * \param iterator iterator for all ghost atoms.
      */
      void begin(ConstGhostIterator& iterator) const;

      //@}
      /// \name Global Atom Counting
      //@{
     
      /**
      * Return maximum number of atoms on all processors.
      *
      * The return value is one greater than the maximum allowed
      * atom id value, i.e. atom ids are assigned values in the
      * range 0 <= id <= totalAtomCapacity-1.
      */
      int totalAtomCapacity() const;

      #ifdef UTIL_MPI
      /**
      * Compute the total number of local atoms on all processors.
      *
      * This is an MPI reduce operation, and thus must be called on
      * all processors. The resulting sum is stored only on the master
      * processor, with rank 0. It may be retrieved by a subsequent 
      * call to nAtomTotal() on the master (rank 0) processor. A null
      * value is stored on all other processors.
      *
      * Upon return, the value is marked as set (i.e., known) on all
      * processors in the communicator. This can be cleared by calling 
      * the unSetNAtomTotal() function on all processors. 
      *
      * If computeNAtomTotal function is called when the value of 
      * nAtomTotal is already set (and thus presumably already known), 
      * the function will return without doing anything.
      *
      * \param communicator MPI communicator for this system.
      */
      void computeNAtomTotal(MPI::Intracomm& communicator);
      #endif

      /**
      * Get total number of atoms on all processors.
      *
      * This function should only be called on the master (rank = 0).
      * The return value is computed by a previous invocation of 
      * computeNAtomTotal(), which must be called on all processors.
      */
      int nAtomTotal() const;

      /**
      * Unset value of NAtomTotal (mark as unknown). 
      *
      * Thus function must be called simultaneously on all processors.
      * It should be called immediately after any operation that changes
      * the number of atoms per processor.
      */
      void unsetNAtomTotal();

      //@}
      /// \name Accessors for Member Objects
      //@{
     
      /**
      * Return the AtomMap by const reference.  
      */
      const AtomMap& map() const;

      #ifdef UTIL_MPI
      /**
      * Get the AtomDistributor by reference.
      */
      AtomDistributor& distributor();

      /**
      * Get the AtomCollector by reference.
      */
      AtomCollector& collector();
      #endif
   
      //@}
      /// \name Miscellaneous Accessors 
      //@{

      /**
      * Has this object been initialized?
      *
      * An AtomStorage is initialized by calling readParam or setParam.
      */
      bool isAllocated() const;

      /**
      * Return true if the container is valid, or throw an Exception.
      */
      bool isValid() const;

      #ifdef UTIL_MPI
      /**
      * Return true if the container is valid, or throw an Exception.
      *  
      * \param communicator communicator for all domain processors.
      */
      bool isValid(MPI::Intracomm& communicator) const;
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
      * Get the maximum number of local atoms encountered thus far.
      *
      * Call only on master.
      */
      int maxNAtom() const;

      /**
      * Get the maximum number of ghost atoms encountered thus far.
      *
      * Call only on master.
      */
      int maxNGhost() const;

      //@}

   private:

      // Array that holds all available local Atom objects.
      AtomArray  atoms_;

      // Set of pointers to local atoms.
      ArraySet<Atom>  atomSet_;

      // Stack of pointers to unused local Atom objects.
      ArrayStack<Atom>  atomReservoir_;

      // Array that holds all available ghost Atom objects.
      AtomArray  ghosts_;

      // Set of pointers to ghost atoms.
      ArraySet<Atom>  ghostSet_;

      // Stack of pointers to unused ghost Atom objects.
      ArrayStack<Atom>  ghostReservoir_;

      // Map of atomIds to atom pointers.
      AtomMap  map_;

      // Array of stored old positions.
      DArray<Vector>  snapshot_;

      // Pointer to space for a new local Atom
      Atom*  newAtomPtr_;

      // Pointer to space for a new ghost Atom.
      Atom*  newGhostPtr_;

      // Capacity for local atoms on this processor.
      int  atomCapacity_;

      // Capacity for ghost atoms on this processors.
      int  ghostCapacity_;

      // Maximum number of atoms on all processors, maximum id + 1
      int  totalAtomCapacity_;

      /// Maximum number of atoms on this proc since stats cleared.
      int  maxNAtomLocal_; 
   
      /// Maximum number of ghosts on this proc since stats cleared.
      int  maxNGhostLocal_;

      #ifdef UTIL_MPI
      // Total number of local atoms on all processors.
      Setable<int>  nAtomTotal_;

      /// Maximum of maxNAtomLocal_ on all procs (defined on master).
      Setable<int>  maxNAtom_;     

      /// Maximum of maxNGhostLocal_ on all procs (defined on master).
      Setable<int>  maxNGhost_;     

      // Distributor and collector
      AtomDistributor distributor_;
      AtomCollector collector_;
      #endif

      // Is addition or removal of atoms forbidden?
      bool locked_;

      // Has memory been allocated to store atoms?
      bool isAllocated_;

      // Are atomic coordinates Cartesian (true) or generalized (false)?
      bool isCartesian_;

      /*
      * Allocate and initialize all private containers.
      */
      void allocate();

      friend class ::AtomStorageTest;
    
   };

   // Inline member function definitions

   inline int AtomStorage::nAtom() const
   { return atomSet_.size(); }

   inline int AtomStorage::nGhost() const
   { return ghostSet_.size(); }

   #ifdef UTIL_MPI
   inline AtomDistributor& AtomStorage::distributor()
   {  return distributor_; }

   inline AtomCollector& AtomStorage::collector()
   {  return collector_; }
   #endif

   inline int AtomStorage::atomCapacity() const
   { return atomCapacity_; }

   inline int AtomStorage::ghostCapacity() const
   { return ghostCapacity_; }

   inline int AtomStorage::totalAtomCapacity() const
   { return totalAtomCapacity_; }

   inline bool AtomStorage::isCartesian() const
   { return isCartesian_; }

   inline const AtomMap& AtomStorage::map() const
   { return map_; }

   inline bool AtomStorage::isAllocated() const
   {  return isAllocated_; }

   /*
   * On master processor (rank=0), stored value of total number of atoms.
   */
   inline int AtomStorage::nAtomTotal() const
   {
      #ifdef UTIL_MPI
      return nAtomTotal_.value();
      #else
      return atomSet_.size();
      #endif
   }

}
#endif
