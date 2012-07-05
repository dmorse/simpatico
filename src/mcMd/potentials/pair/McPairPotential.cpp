#ifndef MCMD_MC_PAIR_POTENTIAL_CPP
#define MCMD_MC_PAIR_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McPairPotential.h"
#include <mcMd/simulation/System.h> 
#include <mcMd/simulation/Simulation.h> 
#include <util/boundary/Boundary.h> 
#include <mcMd/chemistry/Atom.h> 
#include <mcMd/chemistry/Molecule.h> 

#include <util/global.h> 

#include <fstream>

namespace McMd
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   McPairPotential::McPairPotential(System& system)
    : ParamComposite(),
      SubSystem(system)
   {}
 
   /* 
   * Destructor. 
   */
   McPairPotential::~McPairPotential() 
   {}

   /* 
   * Build the CellList.
   */ 
   void McPairPotential::buildCellList() 
   {

      // Adjust number of cells in each direction to current boundary
      cellList_.makeGrid(boundary(), maxPairCutoff());

      // Clear the cellList
      cellList_.clear();

      // Add all atoms to cellList_ 
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               boundary().shift(atomIter->position());
               cellList_.addAtom(*atomIter);
            }
         }
      }

   }

}
#endif
