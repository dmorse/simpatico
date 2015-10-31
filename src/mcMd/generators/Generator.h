#ifndef MCMD_GENERATOR_H
#define MCMD_GENERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>  // function argument, template
#include <util/boundary/Boundary.h>  // typedef

namespace McMd
{

   using namespace Util;

   class Species;
   class Simulation;
   class System;
   class Atom;
   class Molecule;
   class CellList;
   #ifdef INTER_BOND
   class BondPotential;
   #endif

   /**
   * Generates initial configurations for all molecules in a species.
   */
   class Generator 
   {
   public:

      /**
      * Constructor.
      *
      * \param species  associated Species object
      * \param system  parent System
      */
      Generator(Species& species, System& system);

      #ifdef INTER_BOND
      /**
      * Create an association with a BondPotential.
      *
      * \param bondPotential BondPotential
      */
      void setBondPotential(BondPotential& bondPotential);
      #endif

      /**
      * Generate nMolecule molecules of the associated Species.
      *
      * \pre CellList must be allocated.
      * \pre Size of diameters array must equal the number of atom types
      * 
      * \param nMolecule desired number of molecules
      * \param diameters array of excluded volume diameters for atomTypes
      * \param cellList CellList
      */
      virtual bool generate(int nMolecule,
                            const DArray<double>& diameters,
                            CellList& cellList);

      // Static member function
   
      /**
      * Check if cell list is allocated, allocate if necessary.
      * 
      * \param system parent System object
      * \param diameters array of excluded volume diameters
      * \param cellList CellList to be allocated
      */
      static void 
      allocateCellList(System& system, 
                       const DArray<double>& diameters, 
                       CellList& cellList);

   protected:

      /**
      * Attempt to place an atom.
      *
      * Checks if the atom position is within the appropriate cutoff
      * distance of any atoms that are already in the cell list. The
      * cutoff distance for atoms of types i and j is the average of
      * the excluded volume diameters for types i and j. If the atom 
      * position satisfies this geometrical constraint, add the atom
      * to the cell list and return true. If it does not, return false.
      * 
      * \param atom  new Atom, with proposed position already set
      * \param diameters  array of excluded volume diameters for types
      * \param cellist  CellList object storing existing atoms
      * \return true on success, false on failure.
      */
      bool attemptPlaceAtom(Atom& atom, 
                            const DArray<double>& diameters,
                            CellList& cellList);

      /**
      * Attempt to place a molecule (pure virtual).
      *
      * \param molecule new molecule, with unknown atomic positions
      * \param diameters  array of excluded volume diameters for types
      * \param cellist  CellList object storing existing atoms
      * \return true on success, false on failure.
      */
      virtual
      bool attemptPlaceMolecule(Molecule& molecule, 
                                const DArray<double>& diameters,
                                CellList& cellList) = 0;

      /**
      * Get the associated Species by reference.
      */
      const Species& species();

      /**
      * Get the associated Simulation by reference.
      */
      Simulation& simulation();

      /**
      * Get the associated System by reference.
      */
      System& system();

      /**
      * Get the associated Boundary by reference.
      */
      const Boundary& boundary() const;

      #ifdef INTER_BOND
      /**
      * Get the associated BondPotential by reference.
      */
      const BondPotential& bondPotential();
      #endif

   private:

      /// Pointer to associated Species
      const Species* speciesPtr_;

      /// Pointer to associated Simulation
      Simulation* simulationPtr_;

      /// Pointer to associated System
      System* systemPtr_;

      /// Pointer to associated Boundary
      Boundary* boundaryPtr_;

      /// Pointer to associated BondPotential
      #ifdef INTER_BOND
      const BondPotential* bondPotentialPtr_;
      #endif

   };

   // Inline functions

   /*
   * Get the associated Species by reference.
   */
   inline const Species& Generator::species()
   {  return *speciesPtr_; }

   /**
   * Get the associated Simulation by reference.
   */
   inline Simulation& Generator::simulation()
   {  return *simulationPtr_; }

   /*
   * Get the associated System by reference.
   */
   inline System& Generator::system()
   {  return *systemPtr_; }

   /*
   * Get the associated Boundary by reference.
   */
   inline const Boundary& Generator::boundary() const
   {  return *boundaryPtr_; }

   #ifdef INTER_BOND
   /*
   * Get the associated BondPotential by reference.
   */
   inline const BondPotential& Generator::bondPotential()
   {  return *bondPotentialPtr_; }
   #endif

}
#endif
