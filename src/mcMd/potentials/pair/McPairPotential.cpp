/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
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
      SystemInterface(system)
   {  setClassName("McPairPotential"); }
 
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
      // Set up a grid of empty cells.
      cellList_.setup(boundary(), maxPairCutoff());

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
