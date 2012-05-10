#ifndef DDMD_ATOM_STORAGE_H
#define DDMD_ATOM_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <ddMd/chemistry/Atom.h>         // member template parameter
#include <ddMd/chemistry/Group.h>        // Used in template methods
#include <util/containers/DArray.h>      // member template
#include <util/containers/ArraySet.h>    // member template
#include <util/containers/ArrayStack.h>  // member template
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   class Atom;
   class AtomIterator;
   class ConstAtomIterator;
   class GhostIterator;
   class ConstGhostIterator;

   /**
   * A container for all the atoms and ghost atoms on this processor.
   *
   * \ingroup DdMd_Storage_Module
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
      * Read parameters, allocate memory and initialize.
      *
      * Parameters (file format):
      *  - atomCapacity      [int]  max number of atoms owned by processor.
      *  - ghostCapacity     [int]  max number of ghosts on this processor.
      *  - totalatomCapacity [int]  max number of atoms on all processors.
      *
      * \param in input parameter stream.
      */
      virtual void readParam(std::istream& in);

      /**
      * Set parameters, allocate memory and initialize.
      *
      * \param atomCapacity      max number of atoms owned by processor.
      * \param ghostCapacity     max number of ghosts on this processor.
      * \param totalAtomCapacity max number of atoms on all processors.
      */
      void initialize(int atomCapacity, int ghostCapacity,
                    int totalAtomCapacity);

      // Atom Mutators

      /**
      * Returns pointer an address available for a new Atom.
      *
      * This method returns the address of an Atom object that can 
      * be used for a new local Atom. It does not add the Atom the
      * the atom set, and so must be followed by a call to
      * addNewAtom() to do this. 
      *
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
      *
      * \return address for new atom.
      */
      Atom* newAtomPtr(); 

      /**
      * Finalize addition of the most recent new atom.
      *
      * This method adds the atom that was returned by the most 
      * recent call to newAtomPtr to the atom set. The global atom 
      * id must be set before calling this function, by calling 
      * Atom::setId(int), because the algorithm uses the global
      * id returned by Atom::id(). 
      */
      void addNewAtom(); 

      /**
      * Add atom with specified global id.
      * 
      * This method adds a new atom to the atom set and returns
      * a pointer to the address of the new atom. It is equivalent to
      * the following, in which storage is an instance of AtomStorage
      * and ptr is an Atom pointer:
      * \code
      * ptr = storage.newAtomPtr;
      * ptr->setId(id);
      * storage.addNewAtom();
      * \endcode
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

      // Ghost Atom Mutators

      /**
      * Returns pointer an address available for a new ghost Atom.
      *
      * This method returns the address of an Atom object that can 
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
      * This method adds the atom that was returned by the most recent 
      * call to newGhostPtr to the ghost atom set. The global atom 
      * id must be set before calling this function, by calling 
      * Atom::setId(int), because an data structure uses the global
      * id returned by Atom::id(). 
      */
      void addNewGhost(); 

      /**
      * Add ghost atom with specified global id.
      * 
      * This method adds a new atom to the ghost atom set and returns
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

      // Snapshots.

      /**
      * Record current positions of all local atoms.
      * 
      * This method also locks the storage, prohibiting addition
      * or removal of atoms or ghosts until clearSnapshot is called.
      */
      void makeSnapshot();

      /**
      * Clear previous snapshot.
      *
      * This method removes the lock imposed by a previous call
      * to takeSnapshot(), allowing changes to atom and ghost sets.
      */
      void clearSnapshot();

      /**
      * Return max. squared displacement on this processor since last snapshot.
      *
      * Note: This is a local operation, and returns only the maximum on this
      * processor. 
      *
      * Throws exception if no valid snapshot is available.
      */
      double maxSqDisplacement();

      #ifdef UTIL_MPI

      /**
      * Determine whether an atom exchange and reneighboring is needed.
      *
      * \param  communicator domain communicator
      * \param  skin         pair potential skin for pair list
      * \return true if reneighboring is necessary
      */
      bool needExchange(MPI::Intracomm& communicator, double skin);

      #else

      /**
      * Determine whether an atom exchange and reneighboring is needed.
      *
      * \param skin  pair potential skin for pair list
      * \return true if reneighboring is necessary
      */
      bool needExchange(double skin);

      #endif

      // Accessors (iterators)

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

      // Accessors (miscellaneous)

      /**
      * Return pointer to Atom with specified id.
      *
      * This method returns a pointer to an Atom or ghost Atom with
      * the specified id if it is exists on this processor, or returns
      * a null pointer if the atom does not exist on this processor.
      */
      Atom* find(int atomId) const;  

      /**
      * Count the number of atoms in Group and in this AtomStorage.
      *
      * Preconditions: 
      * 1) All atom ids in the Group<N> must be set to valid values,
      *    in the range 0 <= atomId(i) < totalAtomCapacity.
      * 2) No ghost atoms may exist in this AtomStorage. The 
      *    method throws an exception if it finds a ghost atom.
      *
      * \param group Group<N> object with known atom ids. 
      */ 
      template <int N> 
      int countGroupAtoms(Group<N>& group) const;

      /**
      * Set handles to atoms in a Group<N> object.
      *
      * On entry group is a Group<N> object for which the atom
      * ids have been set for all N atoms in the group, but the
      * pointers may not yet have been set. On exit, both pointers
      * and values of atom ownerId are set are set for all atoms
      * that are found in this AtomStorage. Pointers for any atoms 
      * that are not found on this storage are set to null values.
      *
      * Precondition: All atom ids in the Group must be set to
      * values in the range 0 <= atomId(i) < totalAtomCapacity.
      *
      * \param group Group<N> object with known atom ids. 
      * \return number of atoms found on this processor.
      */ 
      template <int N> 
      int findGroupAtoms(Group<N>& group) const;

      /**
      * Return current number of atoms (excluding ghosts)
      */
      int nAtom() const;

      /**
      * Return current number of ghost atoms.
      */
      int nGhost() const;

      /**
      * Return capacity for atoms (excluding ghosts).
      */
      int atomCapacity() const;

      /**
      * Return capacity for ghost atoms
      */
      int ghostCapacity() const;

      /**
      * Return maximum number of atoms on all processors.
      *
      * Atom ids are labelled from 0, ..., totalAtomCapacity-1
      */
      int totalAtomCapacity() const;

      #ifdef UTIL_MPI
      /**
      * Compute and store the total number of atoms on all processors.
      *
      * This is an MPI reduce operation. The correct result is stored and
      * returned only on the rank 0 processor. On other processors, the
      * method stores a null value of -1.
      *
      * \param  communicator MPI communicator for this system.
      * \return on master node, return total number of atoms.
      */
      void computeNAtomTotal(MPI::Intracomm& communicator);
      #endif
   
      /**
      * On master processor (rank=0), stored value of total number of atoms.
      *
      * This method should only be called on the master node (rank = 0).
      * The return value is computed by a previous call of computeNAtomTotal.
      */
      int nAtomTotal() const;

      /**
      * Has this object been initialized?
      *
      * An AtomStorage is initialized by calling readParam or setParam.
      */
      bool isInitialized() const;

      /**
      * Return true if the container is valid, or throw an Exception.
      */
      bool isValid() const;

   private:

      // Array that holds all available local Atom objects.
      DArray<Atom>     atoms_;

      // Set of pointers to local atoms.
      ArraySet<Atom>   atomSet_;

      // Stack of pointers to unused local Atom objects.
      ArrayStack<Atom> atomReservoir_;

      // Array that holds all available local Atom objects.
      DArray<Atom>     ghosts_;

      // Set of pointers to ghost atoms.
      ArraySet<Atom>   ghostSet_;

      // Stack of pointers to unused ghost Atom objects.
      ArrayStack<Atom> ghostReservoir_;

      // Array of pointers to atoms, indexed by Id.
      // Elements corresponding to absent atoms hold null pointers.
      DArray<Atom*>    atomPtrs_;

      // Array of stored old positions.
      DArray<Vector>   snapshot_;

      // Pointer to space for a new local Atom
      Atom* newAtomPtr_;

      // Pointer to space for a new ghost Atom.
      Atom* newGhostPtr_;

      // Capacity for local atoms on this processor.
      int atomCapacity_;

      // Capacity for ghost atoms on this processors.
      int ghostCapacity_;

      // Maximum number of atoms on all processors, maximum id + 1
      int totalAtomCapacity_;

      // Total number of local atoms on all processors.
      int nAtomTotal_;

      // Is addition or removal of atoms forbidden?
      bool locked_;

      // Is this object initialized (has memory been allocated?).
      bool isInitialized_;

      /*
      * Allocate and initialize all private containers.
      */
      void allocate();
    
   };

   // Inline method definitions

   inline int AtomStorage::nAtom() const
   { return atomSet_.size(); }

   inline int AtomStorage::nGhost() const
   { return ghostSet_.size(); }

   inline int AtomStorage::atomCapacity() const
   { return atomCapacity_; }

   inline int AtomStorage::ghostCapacity() const
   { return ghostCapacity_; }

   inline int AtomStorage::totalAtomCapacity() const
   { return totalAtomCapacity_; }

   /*
   * On master processor (rank=0), stored value of total number of atoms.
   */
   inline int AtomStorage::nAtomTotal() const
   {
      #ifdef UTIL_MPI
      if (nAtomTotal_ < 0) UTIL_THROW("Value not set: nAtomTotal < 0");
      return nAtomTotal_;
      #else
      return atomSet_.size();
      #endif
   }

   // Template method definition

   /*
   * Count the number of atoms in this Group and this AtomStorage.
   */
   template <int N>
   int AtomStorage::countGroupAtoms(Group<N>& group) const
   {
      Atom* ptr;
      int nAtom = 0;
      for (int i = 0; i < N; ++i) {
         ptr = atomPtrs_[group.atomId(i)];
         if (ptr) {
            if (ptr->isGhost()) {
               UTIL_THROW("Found ghost atom");
            }
            ++nAtom;
         }
      }
      return nAtom;
   }

   /*
   * Set pointers to atoms in a Group<N> object.
   */
   template <int N>
   int AtomStorage::findGroupAtoms(Group<N>& group) const
   {
      Atom* ptr;
      int nAtom = 0;
      for (int i = 0; i < N; ++i) {
         ptr = atomPtrs_[group.atomId(i)];
         if (ptr) {
            group.setAtomPtr(i, ptr);
            ++nAtom;
         } else {
            group.clearAtomPtr(i);
         }
      }
      return nAtom;
   }

}
#endif
