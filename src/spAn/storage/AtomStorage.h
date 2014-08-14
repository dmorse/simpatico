#ifndef SPAN_ATOM_STORAGE_H
#define SPAN_ATOM_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/chemistry/Atom.h>              // member (template argument)
#include <spAn/chemistry/Group.h>             // member (template argument)

#include <util/containers/DArray.h>           // member (template)
#include <util/containers/DSArray.h>          // member (template)
#include <util/containers/ArrayIterator.h>    // inline function


namespace SpAn 
{

   using namespace Util;

   /**
   * Container for a set of atoms.
   *
   * \ingroup SpAn_Storage_Module
   */
   class AtomStorage 
   {

   public:

      typedef ArrayIterator<Atom> Iterator;

      /**
      * Constructor
      */
      AtomStorage();

      /**
      * Destructor
      */
      ~AtomStorage();

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
      Atom* newPtr();

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
      Atom* ptr(int id);

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
      DSArray<Atom> atoms_;

      /// Pointers to atoms indexed by ids. Missing atoms are null pointers.
      DArray<Atom*> atomPtrs_;

      /// Pointer to new atom.
      Atom* newPtr_;

   };

   // inline functions

   /*
   * Return number of atoms.
   */
   inline int AtomStorage::size() const
   {  return atoms_.size(); }

   /*
   * Get atom capacity (maximum id + 1).
   */ 
   inline
   int AtomStorage::capacity() const
   {  return atoms_.capacity(); }

   /*
   * Return a pointer to an atom with a specific id.
   */
   inline Atom* AtomStorage::ptr(int id)
   {  return atomPtrs_[id]; }

   /*
   * Initialize an iterator for atoms.
   */
   inline 
   void AtomStorage::begin(AtomStorage::Iterator& iter)
   {  atoms_.begin(iter); }

}
#endif
