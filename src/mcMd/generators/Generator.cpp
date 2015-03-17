/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "Generator.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <mcMd/species/Species.h>
#include <mcMd/neighbor/CellList.h>

namespace McMd
{

   using namespace Util;

   Generator::Generator(Species& species, System& system)
    : speciesPtr_(&species),
      simulationPtr_(&system.simulation()),
      systemPtr_(&system),
      boundaryPtr_(&system.boundary()),
      bondPotentialPtr_(0)
   {}

   #ifdef INTER_BOND
   void Generator::setBondPotential(BondPotential& bondPotential)
   {  bondPotentialPtr_ = &bondPotential; }
   #endif

   /*
   * Attempt to place an atom.
   */
   bool Generator::attemptPlaceAtom(Atom& atom,
                             const DArray<double>& diameters, 
                             CellList& cellList)
   {
      boundary().shift(atom.position());
      Vector pos = atom.position();
      double rSq, di, dj, dSq;
      CellList::NeighborArray neighbors;
      Atom* neighborPtr;
      int n;

      cellList.getNeighbors(pos, neighbors);
      di = diameters[atom.typeId()];
      n = neighbors.size();
      for (int j = 0; j < n; ++j) {
         neighborPtr = neighbors[j];
         rSq = boundary().distanceSq(neighborPtr->position(), pos);
         dj = diameters[neighborPtr->typeId()];
         dSq = 0.5*(di + dj);
         dSq *= dSq;
         if (rSq < dSq) {
            return false;
         }
      } 
      cellList.addAtom(atom);
      return true;
   }

   /*
   * Generate random molecules
   */
   bool Generator::generate(int nMolecule, 
                            const DArray<double>& diameters, 
                            CellList& cellList)
   {
      UTIL_CHECK(nMolecule <= species().capacity());
      UTIL_CHECK(cellList.isAllocated());

      // If cell list is not allocated, then allocate.
      if (!cellList.isAllocated()) {
         double maxDiameter = 0.0;
         for (int iType = 0; iType < simulation().nAtomType(); iType++) {
            if (diameters[iType] > maxDiameter)
               maxDiameter = diameters[iType];
         }
         cellList.allocate(simulation().atomCapacity(), 
                           boundary(), maxDiameter);
      }

      // Attempt to place all molecules in Species
      bool success;
      int iMol, iAttempt, maxAttempt;
      maxAttempt = 100;
      Simulation& sim = simulation();
      for (iMol = 0; iMol < nMolecule; ++iMol) {
         Molecule &newMolecule= sim.getMolecule(species().id());
         system().addMolecule(newMolecule);
         success = false;
         iAttempt = 0;
         while (!success && iAttempt < maxAttempt) {
            success = attemptPlaceMolecule(newMolecule, 
                                           diameters, cellList);
            ++iAttempt;
         }
         if (!success) {
            return false;
         }
      }
      return true;
   }

   /*
   * If cell list is not allocated, allocate it.
   */
   void Generator::allocateCellList(System& system,
                                    const DArray<double>& diameters, 
                                    CellList& cellList)
   {
      Simulation& simulation = system.simulation();
      UTIL_CHECK(diameters.capacity() == simulation.nAtomType());
      if (!cellList.isAllocated()) {
         double maxDiameter = 0.0;
         for (int iType = 0; iType < diameters.capacity(); iType++) {
            if (diameters[iType] > maxDiameter) {
               maxDiameter = diameters[iType];
            }
         }
         cellList.allocate(simulation.atomCapacity(), 
                           system.boundary(), maxDiameter);
      }
   }

}
