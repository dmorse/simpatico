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
      * On entry, group is a Group<N> object for which the atom ids
      * for all N atoms in the Group have been set to valid values, 
      * in the range 0 <= atomId < totalAtomCapacity, but in which
      * some or all pointers have not been set. 
      *
      * On exit, pointers are set for all atoms are present in this 
      * AtomMap, or set to null for absent atoms. This method overwrites
      * all old pointer values.
      *
      * \param group Group<N> object with known atom ids. 
      * \return number of atoms found on this processor.
      */ 
      template <int N> 
      int findGroupAtoms(Group<N>& group) const;

      /**
      * Set handles to atoms in a Group<N> object.
      *
      * On entry, group is a Group<N> object for which the atom ids
      * for all N atoms in the Group have been set to valid values, 
      * in the range 0 <= atomId < totalAtomCapacity, but in which
      * some or all pointers have not been set. The AtomMap may not
      * contain any ghosts.
      *
      * On exit, pointers are set for all local atoms that exist in 
      * this AtomMap, or set to null otherwise. All old pointer values 
      * are overwritten. 
      *
      * \param group Group<N> object with known atom ids. 
      * \return number of atoms found on this processor.
      */ 
      template <int N> 
      int findGroupLocalAtoms(Group<N>& group) const;

      /**
      * Set handles to ghost atoms in a Group<N> object.
      *
      * On entry, group is a Group<N> object for which the atom ids
      * for all N atoms in the Group have been set to valid values, 
      * in the range 0 <= atomId < totalAtomCapacity, and in which
      * all pointers to local atoms have been set, but in which some
      * or all pointers to ghosts have not been set. This function
      * should be called after this AtomMap contains all ghost atoms.
      * The function does not modify non-null pointers set previously.
      *
      * On exit, pointers are set for all ghost atoms present in this
      * AtomMap.
      *
      * \param group Group<N> object with known atom ids. 
      * \return number of atoms found on this processor.
      */ 
      template <int N> 
      int findGroupGhostAtoms(Group<N>& group) const;

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
   {  return nLocal_; }

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

   /*
   * Set pointers to all atoms in a Group<N> object.
   */
   template <int N>
   int AtomMap::findGroupLocalAtoms(Group<N>& group) const
   {
      Atom* ptr;
      int nAtom = 0;
      for (int i = 0; i < N; ++i) {
         ptr = atomPtrs_[group.atomId(i)];
         if (ptr) {
            assert(!ptr->isGhost());
            assert(ptr->atomId() == group.atomId(i));
            group.setAtomPtr(i, ptr);
            ++nAtom;
         } else {
            group.clearAtomPtr(i);
         }
      }
      return nAtom;
   }

   /*
   * Set pointers to atoms in a Group<N> object.
   */
   template <int N>
   int AtomMap::findGroupGhostAtoms(Group<N>& group) const
   {
      Atom* ptr;
      int nAtom = 0;
      for (int i = 0; i < N; ++i) {
         if (group.atomPtr(i)) {
            ++nAtom;
         } else {
            ptr = atomPtrs_[group.atomId(i)];
            if (ptr) {
               assert(ptr->isGhost());
               assert(ptr->atomId() == group.atomId(i));
               group.setAtomPtr(i, ptr);
               ++nAtom;
            }
         }
      }
      return nAtom;
   }

}
#endif
