#ifndef DDMD_MASK_H
#define DDMD_MASK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace DdMd
{

   class Atom;

   /**
   * Set of Atoms for which pair interactions with a target Atom are "masked".
   *
   * A Mask stores int identifiers for a set of Atoms for which the non-bonded
   * pair interactions with a single target atom are suppressed, or "masked".
   * Each Mask object is associated with one target Atom. 
   *
   * A Mask could, for example, contain all atoms that are directly connected 
   * to the target atom by 2-body covalent bonds.
   *
   * \ingroup DdMd_Chemistry_Module
   */
   class Mask 
   {

   public:

      /// Maximum number of masked atoms per target atom.
      static const int Capacity  = 4;

      /**
      * Constructor.
      */
      Mask();

      /**
      * Clear the masked set (remove all atoms).
      */
      void clear();

      /**
      * Add an Atom to the masked set.
      *  
      * \param id global index (tag) of atom to be added
      */
      void append(int id);

      /**
      * True if the atom is in the masked set for the target Atom.
      *
      * \param  id integer id of atom to be tested
      * \return true if atom is masked, false otherwise
      */
      bool isMasked(int id) const;

      /**
      * Return value of atom index number i.
      *
      * \param i array index for desired atom index
      */
      int operator [] (int i) const;

      /**
      * Return the number of masked atoms.
      */
      int size() const;

   private:

      /// Integer ids to of masked Atoms.
      int atomIds_[Capacity];      

      /// Number of masked Atoms.
      int  size_;                 

   }; 

   // Inline methods

   // Check if an Atom is in the masked set.
   inline bool Mask::isMasked(int id) const
   {
      for (int i=0; i < size_ ; ++i) {
         if (atomIds_[i] == id) return true;
      }
      return false;
   }

   /*
   * Return the number of masked atoms.
   */
   inline int Mask::size() const
   { return size_; }

   /*
   * Return value of index
   */
   inline int Mask::operator[] (int i) const
   {
      assert(i >= 0);
      assert(i <  size_);
      return atomIds_[i]; 
   }

} 
#endif
