#ifndef DDMD_ATOM_MAP_H
#define DDMD_ATOM_MAP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/chemistry/Atom.h>     // member template argument
#include <util/containers/DArray.h>  // member template
#include <ddMd/chemistry/Group.h>    // member function template
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Associative container for finding atoms identified by integer id.
   *
   * \ingroup DdMd_Storage_Module
   */
   class AtomMap 
   {

   public:

      /**
      * Constructor.
      */
      AtomMap();

      /**
      * Destructor.
      */
      ~AtomMap();

      /**
      * Set parameters, allocate memory and initialize.
      *
      * Call this or (read|load)Parameters to initialize, but not both.
      *
      * \param totalAtomCapacity max number of atoms on all processors.
      */
      void allocate(int totalAtomCapacity);

      /**
      * Add local atom.
      * 
      * \param ptr Pointer to new Atom.
      */
      void addLocal(Atom* ptr); 

      /**
      * Remove a specific Atom.
      *
      * \throw Exception if atom is not present.
      *
      * \param ptr Pointer to Atom to be removed.
      */
      void removeLocal(Atom* ptr); 

      /**
      * Add ghost atom.
      * 
      * \param ptr Pointer to new Atom.
      */
      void addGhost(Atom* ptr);

      /**
      * Remove a ghost Atom.
      *
      * This method throws an exception if no atom with this
      * id is present, but not if it does not match this pointer.
      *
      * \param ptr Pointer to Atom to be removed.
      */
      void removeGhost(Atom* ptr); 

      //@}
      /// \name Accessors 
      //@{

      /**
      * Return pointer to Atom with specified id.
      *
      * This method returns a pointer to an Atom with the specified
      * id if it is present, or returns a null pointer otherwise.
      *
      * \param atomId integer index of atom
      */
      Atom* find(int atomId) const;  

      /**
      * Return the number of local atoms.
      */ 
      int nLocal() const;

      /**
      * Return the number of ghosts with distinct ids.
      */ 
      int nGhost() const;

      /**
      * Set handles to atoms in a Group<N> object.
      *
      * On entry group is a Group<N> object for which the atom
      * ids have been set for all N atoms in the group, but the
      * pointers may not yet have been set. On exit, both pointers
      * and values of atom ownerId are set are set for all atoms
      * that are found in this AtomMap. Pointers for any atoms 
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
      * Check validity of this AtomMap.
      *
      * Returns true if all is ok, or throws an Exception.
      */
      bool isValid() const;

      //@}

   private:

      // Array of pointers to atoms, indexed by Id.
      // Elements corresponding to absent atoms hold null pointers.
      DArray<Atom*>  atomPtrs_;

      /// Number of local atoms in this map.
      int nLocal_;

      /// Number of ghost atoms in this map.
      int nGhost_;

      // Maximum number of atoms on all processors, maximum id + 1
      int totalAtomCapacity_;

      // Has this map been initialized (i.e., allocated)?
      bool isInitialized_;

   };

   // Inline method definitions

   /*
   * Return pointer to Atom with specified id.
   */
   inline Atom* AtomMap::find(int atomId) const
   {  return atomPtrs_[atomId]; }

   /*
   * Return the number of local atoms.
   */ 
   inline int AtomMap::nLocal() const
   { return nLocal_; }

   /*
   * Return the number of ghosts with distinct ids.
   */ 
   inline int AtomMap::nGhost() const
   { return nGhost_; }

   // Template method definition

   /*
   * Set pointers to atoms in a Group<N> object.
   */
   template <int N>
   int AtomMap::findGroupAtoms(Group<N>& group) const
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
