#ifndef MCMD_GENERATOR_H
#define MCMD_GENERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/neighbor/CellList.h>  // member
#include <util/containers/DArray.h>  // function argument, template
#include <util/boundary/Boundary.h>  // typedef

namespace McMd
{

   using namespace Util;

   class Species;
   class Simulation;
   class System;
   #ifdef INTER_BOND
   class BondPotential;
   #endif

   class Generator {

   public:

      Generator(Species& species, System& system);

      #ifdef INTER_BOND
      void setBondPotential(BondPotential& bondPotential);
      #endif

      virtual bool generate(int nMolecule,
                            const DArray<double> diameters,
                            CellList& cellList);

   protected:

      /**
      * Attempt to place an atom.
      *
      * If the atom satisfies the geometrical constraints, add the
      * atom to the cell list and return true. Otherwise return false.
      */
      bool attemptPlaceAtom(Atom& atom, 
                            const DArray<double> diameters,
                            CellList& cellList);

      /**
      * Attempt place a molecule.
      *
      * \return true on success, false on failure.
      */
      virtual
      bool attemptPlaceMolecule(Molecule& molecule, 
                                const DArray<double> diameters,
                                CellList& cellList) = 0;

      const Species& species()
      {  return *speciesPtr_; }

      Simulation& simulation()
      {  return *simulationPtr_; }

      System& system()
      {  return *systemPtr_; }

      const Boundary& boundary() const
      {  return *boundaryPtr_; }

      #ifdef INTER_BOND
      BondPotential& bondPotential()
      {  return *bondPotentialPtr_; }
      #endif

   private:

      const Species* speciesPtr_;

      Simulation* simulationPtr_;

      System* systemPtr_;

      Boundary* boundaryPtr_;

      #ifdef INTER_BOND
      BondPotential* bondPotentialPtr_;
      #endif

   };

}
#endif
