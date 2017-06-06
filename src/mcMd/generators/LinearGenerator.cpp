/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinearGenerator.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/neighbor/CellList.h>
#include <simp/species/Species.h>
#include <util/boundary/Boundary.h>

namespace McMd
{

   class CellList;

   using namespace Util;
   using namespace Simp;

   LinearGenerator::LinearGenerator(Species& species, System& system)
    : Generator(species, system)
   {}

   /*
   * Recursive function to try to place an atom.
   */
   bool 
   LinearGenerator::attemptPlaceMolecule(Molecule& molecule,
                                        Array<double> const & diameters,
                                        CellList& cellList)
   {
      Random& random = simulation().random();

      // Attempt to place first atom at random
      Atom* atomPtr = &molecule.atom(0);
      int maxAttempt = 500;
      int iAttempt = 0;
      bool success = false;
      while (!success && iAttempt < maxAttempt) {
         boundary().randomPosition(random, atomPtr->position());
         success = attemptPlaceAtom(*atomPtr, diameters, cellList);
      }
      if (!success) {
         return false; 
      }

      Vector v;
      //double temperature = system().energyEnsemble().temperature();
      //double beta = 1.0/temperature;
      double beta = 1.0;
      Atom* prevPtr = 0;
      int bondType = 0;
      maxAttempt = 100;
      int nAtom = species().nAtom();
      for (int iAtom = 1; iAtom < nAtom; ++iAtom) {
         atomPtr = &molecule.atom(iAtom);
         prevPtr = &molecule.atom(iAtom-1);
         iAttempt = 0;
         success = false;
         while (!success && iAttempt < maxAttempt) {
            atomPtr->position() = prevPtr->position();
            random.unitVector(v);
            v *= bondPotential().randomBondLength(&random, beta, 
                                                  bondType);
            atomPtr->position() += v;
            success = attemptPlaceAtom(*atomPtr, diameters, cellList);
            ++iAttempt;
         }
         if (!success) {
            for (int j = 0; j < iAtom; ++j) {
               cellList.deleteAtom(molecule.atom(j));
            }
            return false; // Failure to insert an atom
         }
      }

      // Normal termination implies successful molecule insertion.
      return true;
   }

}
