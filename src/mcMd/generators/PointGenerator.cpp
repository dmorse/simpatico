/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PointGenerator.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
//#include <mcMd/neighbor/CellList.h>
#include <util/boundary/Boundary.h>

namespace McMd
{

   class CellList;
   using namespace Util;

   PointGenerator::PointGenerator(Species& species, System& system)
    : Generator(species, system)
   {}

   /*
   * Recursive function to try to place an atom.
   */
   bool 
   PointGenerator::attemptPlaceMolecule(Molecule& molecule,
                                        const DArray<double>& diameters,
                                        CellList& cellList)
   {
      Atom& atom = molecule.atom(0);
      bool success = false;
      int iAttempt = 0;
      int maxAttempt = 1;
      while (!success && iAttempt < maxAttempt) {
         boundary().randomPosition(simulation().random(), atom.position());
         success = attemptPlaceAtom(atom, diameters, cellList);
      }
      return success;
   }

}
