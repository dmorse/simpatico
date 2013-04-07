#ifndef DDMD_MASK_H
#define DDMD_MASK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace DdMd
{

   class Atom;

   /**
   * Set of Atoms for which pair interactions with a parent Atom are "masked".
   *
   * Each Mask object is associated with one parent Atom. The Mask stores 
   * the integer identifiers for a set of Atoms for which non-bonded pair 
   * interactions with the parent atom are suppressed, or "masked". These are
   * generally atoms that are directly bonded to the parent Atom, or part of
   * the same Angle or Dihedral group. The Mask is used during construction 
   * of a Velet pair list to identify nearby atoms for which pair interactions 
   * are suppressed.
   *
   * \ingroup DdMd_Chemistry_Module
   */
   class Mask 
   {

   public:

      /*
      * Maximum number of masked atoms per parent atom.
      */
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
      * True if the atom is in the masked set for the parent Atom.
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

   /*
   * Check if an Atom is masked.
   */
   inline bool Mask::isMasked(int id) const
   {
      for (int i=0; i < size_ ; ++i) {
         if (atomIds_[i] == id) return true;
      }
      return false;
   }

   /*
   * Return the number of masked atoms for this parent Atom.
   */
   inline int Mask::size() const
   { return size_; }

   /*
   * Return value of atom index.
   */
   inline int Mask::operator[] (int i) const
   {
      assert(i >= 0);
      assert(i <  size_);
      return atomIds_[i]; 
   }

} 
#endif
