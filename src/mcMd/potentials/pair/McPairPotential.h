#ifndef MCMD_MC_PAIR_POTENTIAL_H
#define MCMD_MC_PAIR_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>             // base class
#include <mcMd/simulation/SystemInterface.h>       // base class
#include <mcMd/potentials/pair/PairPotential.h>    // base class
#include <mcMd/neighbor/CellList.h>                // member

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   class System;

   /**
   * A PairPotential for MC simulations (abstract).
   *
   * \ingroup McMd_Pair_Module
   */
   class McPairPotential : public ParamComposite, public PairPotential, 
                           public SystemInterface
   {

   public:

      /**   
      * Constructor.
      */
      McPairPotential(System& system);

      /** 
      * Destructor.
      */
      virtual ~McPairPotential();

      /// \name Energy, Force, Stress evaluators (pure virtual)
      //@{

      /**
      * Calculate the nonbonded pair energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return nonbonded pair potential energy of atom
      */
      virtual double atomEnergy(const Atom& atom) const = 0;

      /**
      * Calculate the nonbonded pair energy for an entire Molecule.
      *
      * The return value is the change in total pair potential upon
      * inserting or removing a molecule. This is the sum of all
      * inter-molecular pair interactions involving atoms of this
      * molecule, plus an intramolecular contribution in which each 
      * pair is counted only once.
      *
      * \param  molecule Molecule object of interest
      * \return nonbonded pair potential energy of molecule
      */
      virtual double moleculeEnergy(const Molecule& molecule) const = 0;

      //@}
      /// \name Cell List Management
      //@{

      /**
      * Build the CellList with current configuration.
      *
      * Calls CellList::clear() to clear the CellList,
      * then adds every Atom in this System. Each Atom
      * position is shifted into the primary box by
      * Boundary::shift() before being added.
      */
      void buildCellList();

      /**
      * Add an Atom to the CellList.
      *
      * \param atom Atom object of interest.
      */
      void addAtom(Atom& atom);

      /**
      * Remove an Atom from the CellList.
      *
      * \param atom Atom object of interest
      */
      void deleteAtom(Atom& atom);

      /**
      * Update the cell list to reflect a new position.
      *
      * Use this to update the cell list after an Atom position
      * has already been assigned a new value.
      *
      * \param atom Atom object whose position has been modified.
      */
      void updateAtomCell(Atom &atom);

      /**
      * Move an Atom position, and update the CellList.
      *
      * On output, atom.position() is reset to to position, and the
      * new position is recorded in the CellList. This is equivalent 
      * to atom.position() = position and then updateAtomCell(atom).
      *
      * \param atom     Atom object of interest
      * \param position new position vector
      */
      void moveAtom(Atom& atom, const Vector &position);

      /** 
      * Get the cellList by const reference.
      */
      const CellList& cellList() const;

      //@}

   protected:

      /// Array to hold neighbors returned by a CellList.
      mutable CellList::NeighborArray neighbors_;

      /// Cell list for atom positions.
      CellList cellList_;

   };

   // Inline functions
  
   // Add an atom to CellList.
   inline void McPairPotential::addAtom(Atom &atom)
   {  cellList_.addAtom(atom); }

   // Delete an atom from the CellList.
   inline void McPairPotential::deleteAtom(Atom &atom)
   {  cellList_.deleteAtom(atom); }

   // Update the cell list to reflect a new Atom position.
   inline void McPairPotential::updateAtomCell(Atom &atom)
   {  cellList_.updateAtomCell(atom, atom.position()); }

   // Move atom to a new position.
   inline void McPairPotential::moveAtom(Atom &atom, const Vector &position)
   {
      atom.position() = position;
      cellList_.updateAtomCell(atom, position);
   }

   // Get the cellList by const reference.
   inline const CellList& McPairPotential::cellList() const
   { return cellList_; }

} 
#endif
