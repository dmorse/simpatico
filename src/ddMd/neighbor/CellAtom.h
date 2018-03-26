#ifndef DDMD_CELL_ATOM_H
#define DDMD_CELL_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace DdMd
{

   /**
   * Data for an atom stored within a CellList.
   *
   * A CellList has an array of CellAtom objects, in which atoms in 
   * the same cell are stored in consecutive elements of the array. 
   * A CellAtom has a pointer to the associated Atom object, and also
   * stores a local copy of the atom position and id. The position 
   * and id are set by calling the update() method, which copies the
   * values from values stored in the associated Atom object.
   *
   * Design rationale: In code that uses a CellList to construct a
   * Verlet pair list, the code that find neighbors of a primary atom
   * requires access to the Mask associated with the primary atom, and
   * requires access to the position and atom id of each secondary atom.
   * Storing the position and id of the atom in the CellAtom makes it
   * possible to write the inner loop over secondary atoms as a loop
   * over CellAtom objects using consecutive access, without needing
   * to dereference the pointer to the Atom object. 
   *
   * Usage: The setPtr function should be called in the code that
   * builds a CellList, when atoms are assigned to a cell. At this
   * point the position may be expressed in generalized coordinates.
   * The update() function should later be called for all atoms in 
   * the cell list by the code that builds the pair list, after the
   * atomic positions have been transformed back to Cartesian
   * coordinates (to simplify computation of inter-atomic distances),
   * and just prior to the loop that constructs a pair list. 
   */
   class CellAtom
   {

   public:

      /**
      * Set the pointer to the associated atom.
      *
      * \param atomPtr pointer to Atom.
      */
      void setPtr(Atom* atomPtr);

      /**
      * Update local copies of the atom position and id.
      */
      void update();

      /**
      * Return a pointer to the associated atom.
      */
      Atom* ptr() const;

      /**
      * Return a pointer to the Mask of the associated atom.
      */
      Mask* maskPtr() const;

      /**
      * Return the stored atom position.
      */
      const Vector& position() const;

      /**
      * Return the stored atom id.
      */
      int id() const;

   private:

      // Local copy of atomic position.
      Vector position_;

      // Pointer to associated Atom object.
      Atom* ptr_;

      // Local copy of atom id.
      int id_;

   };

   // Inline functions

   // Set the pointer to the associated atom.
   inline void CellAtom::setPtr(Atom* atomPtr) 
   {  ptr_ = atomPtr; }

   // Update stored copies of the position and id.
   inline void CellAtom::update() 
   {
      position_ = ptr_->position();
      id_ = ptr_->id();
   }

   // Return a pointer to the associated atom.
   inline Atom* CellAtom::ptr() const
   {  return ptr_; }

   // Return a pointer to the Mask of the associated atom.
   inline Mask* CellAtom::maskPtr() const
   {  return &(ptr_->mask()); }

   // Return the stored atom position.
   inline const Vector& CellAtom::position() const
   {  return position_; }

   // Return the stored atom id.
   inline int CellAtom::id() const
   {  return id_; }

}
#endif
