#ifndef MCMD_MASK_H
#define MCMD_MASK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   class Atom;

   /**
   * Set of Atoms for which pair interactions with a target Atom are "masked".
   *
   * An Mask stores identifiers for a set of Atoms for which the non-bonded
   * pair interactions with a single target atom are suppressed, or "masked".
   * Each Mask object is associated with one target Atom. 
   *
   * A Mask could, for example, contain all atoms that are directly connected 
   * to the target atom by 2-body covalent bonds.
   *
   * \ingroup McMd_Chemistry_Module
   */
   class Mask 
   {

   public:

      /**
      * Constructor.
      */
      Mask();

      /**
      * Clear the mask set (remove all atoms).
      */
      void clear();

      /**
      * Add an Atom to the masked set.
      *  
      * \param atom Atom to be added to the masked set.
      */
      void append(const Atom& atom);

      /**
      * True if the atom is in the masked set for the target Atom.
      *
      * \param  atom Atom object to be tested.
      * \return true if atom is masked, false otherwise.
      */
      bool isMasked(const Atom& atom) const;

      /**
      * Return the number of masked atoms.
      */
      int size() const;

   private:

      /// Maximum number of masked atoms per target atom.
      static const int Capacity  = 4;

      /// Pointers to of masked Atoms.
      const Atom* atomPtrs_[Capacity];      

      /// Number of masked Atoms.
      int  size_;                 

      /// Copy constructor. Private to prevent copying.
      Mask(const Mask& other);

   }; 

   // Inline member functions.

   /*
   * Check if an Atom is in the masked set.
   */
   inline bool Mask::isMasked(const Atom& atom) const
   {
      const Atom* ptr = &atom;
      for (int i=0; i < size_ ; ++i) {
         if (atomPtrs_[i] == ptr) return true;
      }
      return false;
   }

   /*
   * Return the number of masked atoms.
   */
   inline int Mask::size() const
   { return size_; }

} 
#endif
