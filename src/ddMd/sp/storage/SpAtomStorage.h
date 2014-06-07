#ifndef DDMD_SP_ATOM_STORAGE_H
#define DDMD_SP_ATOM_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/sp/chemistry/SpAtom.h>              // member (template argument)
#include <ddMd/sp/chemistry/SpGroup.h>             // member (template argument)

#include <util/containers/DArray.h>           // member (template)
#include <util/containers/DSArray.h>          // member (template)
#include <util/containers/ArrayIterator.h>    // inline function


namespace DdMd 
{

   using namespace Util;

   /**
   * Container for a set of atoms.
   *
   * \ingroup DdMd_Sp_Storage_Module
   */
   class SpAtomStorage 
   {

   public:

      typedef ArrayIterator<SpAtom> Iterator;

      /**
      * Constructor
      */
      SpAtomStorage();

      /**
      * Destructor
      */
      ~SpAtomStorage();

      /**
      * Allocate and initialize memory.
      *
      * \param capacity Maximum number of atoms
      */
      void allocate(int capacity);

      /**
      * Return pointer to location for new atom.
      *
      * \return pointer to location of new atom
      */
      SpAtom* newPtr();

      /**
      * Finalize addition of atom (allows lookup by id).
      */
      void add();

      /**
      * Clear all atoms and groups.
      */
      void clear();
  
      /**
      * Get a pointer to an atom by global id.
      */
      SpAtom* ptr(int id);

      /**
      * Initialize an iterator for atoms.
      */
      void begin(Iterator& iter);

      /**
      * Get atom capacity (maximum id + 1).
      */ 
      int capacity() const;

      /**
      * Get number of atoms.
      */ 
      int size() const;

   private:
     
      /// Array of atom objects, added in order read from file.
      DSArray<SpAtom> atoms_;

      /// Pointers to atoms indexed by ids. Missing atoms are null pointers.
      DArray<SpAtom*> atomPtrs_;

      /// Pointer to new atom.
      SpAtom* newPtr_;

   };

   // inline functions

   /*
   * Return number of atoms.
   */
   inline int SpAtomStorage::size() const
   {  return atoms_.size(); }

   /*
   * Get atom capacity (maximum id + 1).
   */ 
   inline
   int SpAtomStorage::capacity() const
   {  return atoms_.capacity(); }

   /*
   * Return a pointer to an atom with a specific id.
   */
   inline SpAtom* SpAtomStorage::ptr(int id)
   {  return atomPtrs_[id]; }

   /*
   * Initialize an iterator for atoms.
   */
   inline 
   void SpAtomStorage::begin(SpAtomStorage::Iterator& iter)
   {  atoms_.begin(iter); }

}
#endif
