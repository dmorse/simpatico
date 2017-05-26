#ifndef MCMD_GENERATOR_H
#define MCMD_GENERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>   // function argument, template
#include <util/boundary/Boundary.h>  // typedef

namespace Simp {
   class Species;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   class Simulation;
   class System;
   class Atom;
   class Molecule;
   class CellList;
   #ifdef SIMP_BOND
   class BondPotential;
   #endif

   /**
   * Generates initial configurations for molecules of one species.
   *
   * \ingroup McMd_Generator_Module
   */
   class Generator 
   {
   public:

      // Non-static member functions

      /**
      * Constructor.
      *
      * \param species  associated Species object
      * \param system  parent System
      */
      Generator(Species& species, System& system);

      /**
      * Destructor.
      */
      virtual ~Generator();

      #ifdef SIMP_BOND
      /**
      * Create an association with a BondPotential.
      *
      * \param bondPotential  BondPotential
      */
      void setBondPotential(BondPotential& bondPotential);
      #endif

      /**
      * Generate nMolecule molecules of the associated Species.
      *
      * \pre Size of diameters array == number of atom types
      * 
      * \param nMolecule  desired number of molecules
      * \param diameters  array of excluded volume diameters for atomTypes
      * \param cellList  CellList object
      */
      virtual bool generate(int nMolecule,
                            Array<double> const & diameters,
                            CellList& cellList);

      // Static member function

      /**
      * Allocate any required memory for the cell list.
      * 
      * \param atomCapacity  maximum allowed atom id + 1
      * \param boundary  Boundary object, periodic unit cell
      * \param diameters  array of excluded volume diameters
      * \param cellList  CellList to be allocated
      */
      static
      void setupCellList(int atomCapacity, 
                         Boundary& boundary,
                         const Array<double>& diameters, 
                         CellList& cellList);

   protected:

      /**
      * Attempt to place an atom.
      *
      * The proposed position for the atom should be set upon entry.
      * This function checks if the atom position is within the 
      * hard-core cutoff distance of any atoms that are already in 
      * the cell list. The cutoff distance for atoms of types i and 
      * j is the average of the excluded volume diameters for types 
      * i and j. If the atom position satisfies this geometrical 
      * constraint, the atom is added to the cell list and returns 
      * true. If it does not, the function returns false.
      * 
      * \param atom  new Atom, with proposed position already set
      * \param diameters  array of excluded diameters for atom types
      * \param cellList  CellList object containing existing atoms
      * \return true on success, false on failure.
      */
      bool attemptPlaceAtom(Atom& atom, 
                            const Array<double>& diameters,
                            CellList& cellList);

      /**
      * Attempt to insert an entire molecule (pure virtual).
      *
      * \param molecule  new molecule, with unknown atomic positions
      * \param diameters  array of excluded volume diameters for types
      * \param cellList  CellList object storing existing atoms
      * \return true on success, false on failure.
      */
      virtual
      bool attemptPlaceMolecule(Molecule& molecule, 
                                const Array<double>& diameters,
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

      #ifdef SIMP_BOND
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
      #ifdef SIMP_BOND
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

   #ifdef SIMP_BOND
   /*
   * Get the associated BondPotential by reference.
   */
   inline const BondPotential& Generator::bondPotential()
   {  return *bondPotentialPtr_; }
   #endif

}
#endif
