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

#ifdef UTIL_CXX11
#include <unordered_map>
#else
#include <map>
#endif

namespace Util {
   template <typename T> class ArraySet;
}

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
      * This function throws an exception if no atom with this
      * id is present, but not if it does not match this pointer.
      *
      * \param ptr Pointer to Atom to be removed.
      */
      void removeGhost(Atom* ptr); 

      /**
      * Clear all ghosts from this map.  
      *
      * \param ghostSet Set containing all ghosts on this processor.
      */
      void clearGhosts(const ArraySet<Atom>& ghostSet);

      //@}
      /// \name Accessors 
      //@{

      /**
      * Return pointer to Atom with specified id.
      *
      * This function returns a pointer to an Atom with the specified
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
      int nGhostDistinct() const;

      /**
      * Return the number of ghosts, including images.
      */ 
      int nGhost() const;

      /**
      * Set handles to local atoms in a Group<N> object.
      *
      * On entry, group is a Group<N> object for which the atom ids
      * for all N atoms in the Group have been set to valid values, 
      * in the range 0 <= atomId < totalAtomCapacity, but in which
      * some or all pointers have not been set. The AtomMap may not
      * contain any ghosts. 
      *
      * On exit, pointers are set correctly for all local atoms that 
      * exist in this AtomMap, or set to null for absent atoms. All 
      * old pointer values are overwritten. 
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
      * all pointers to local atoms have been set, but in which no
      * pointers to ghosts have been set. This function may only be
      * called after all pointers have been set for all local atoms
      * and after this AtomMap contains all ghost atoms. 
      *
      * On exit, pointers are set for all ghost atoms present in
      * this AtomMap.
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

      #ifdef UTIL_CXX11
      typedef std::unordered_multimap<int, Atom*> GhostMap;
      #else
      typedef std::multimap<int, Atom*> GhostMap;
      #endif

      // Array of pointers to atoms, indexed by Id.
      // Elements corresponding to absent atoms hold null pointers.
      DArray<Atom*> atomPtrs_;

      // Map for extra ghost images
      GhostMap ghostMap_;

      /// Number of local atoms in this map.
      int nLocal_;

      /// Number of ghost atom with distinct ids in this map.
      int nGhostDistinct_;

      // Maximum number of atoms on all processors, maximum id + 1
      int totalAtomCapacity_;

      // Has this map been initialized (i.e., allocated)?
      bool isInitialized_;

      /*
      *  Design / invariants:
      *
      *  - If a local atom with id i is present, atomPtrs_[i]
      *    contains a pointer to that atom.
      *
      *  - If one or more ghosts with an atom id i are present,
      *    but there is no local atom with that id, atomPtrs_[i]
      *    contains a pointer to one such ghost.
      *
      *  - ghostMap_ contains pointers to all ghosts except those
      *    in atomPtrs_, stored using atom indices as keys. 
      *
      * One image of each physical atom, identified by id, is thus 
      * stored * in atomPtrs_, while ghostMap_ holds any "extra"
      * ghost images of atoms. If this processor does not contain
      * multiple images of any particle, ghostMap_ will be empty.
      */

   };

   // Inline function definitions

   /*
   * Return pointer to an Atom with specified id.
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
   inline int AtomMap::nGhostDistinct() const
   { return nGhostDistinct_; }

   /*
   * Return the total number of ghosts, including images.
   */ 
   inline int AtomMap::nGhost() const
   { return nGhostDistinct_ + ghostMap_.size(); }

   // Function template definitions

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
            assert(ptr->id() == group.atomId(i));
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
         if (group.atomPtr(i) != 0) {
            // assert(!(ptr->isGhost()));
            ++nAtom;
         } else {
            int atomId = group.atomId(i);
            ptr = atomPtrs_[atomId];
            if (ptr) {
               assert(ptr->isGhost());
               assert(ptr->id() == atomId);
               group.setAtomPtr(i, ptr);
               ++nAtom;
            } else {
               group.clearAtomPtr(i);
            }
         }
      }
      return nAtom;
   }

}
#endif
