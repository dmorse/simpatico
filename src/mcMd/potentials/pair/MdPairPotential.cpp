#ifndef INTER_NOPAIR
#ifndef MD_PAIR_INTERACTION
#define MD_PAIR_INTERACTION

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdPairPotential.h"
#include <mcMd/simulation/System.h> 
#include <mcMd/simulation/Simulation.h> 
#include <util/boundary/Boundary.h> 

#include <util/global.h> 

#include <fstream>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   MdPairPotential::MdPairPotential(System& system)
    : ParamComposite(),
      SubSystem(system)
   {  setClassName("MdPairPotential"); }
 
   /* 
   * Destructor. 
   */
   MdPairPotential::~MdPairPotential() 
   {}

   /* 
   * Build the PairList.
   */ 
   void MdPairPotential::buildPairList() 
   {

      // Recalculate the grid for the internal CellList
      pairList_.makeGrid(boundary());

      // Clear all atoms from the internal CellList
      pairList_.clear();

      // Add every atom in this System to the CellList
      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            molIter->begin(atomIter); 
            for ( ; atomIter.notEnd(); ++atomIter) {
               #ifdef MCMD_SHIFT
               boundary().shift(atomIter->position(), atomIter->shift());
               #else
               boundary().shift(atomIter->position());
               #endif
               pairList_.addAtom(*atomIter);
            }
         }
      }

      // Use the completed CellList to build the PairList 
      pairList_.build(boundary());
   }

   /* 
   * Clear the PairList statistical accumulators
   */ 
   void MdPairPotential::clearPairListStatistics() 
   {  pairList_.clearStatistics(); }

   /* 
   * Return true if the pair list is current, false if it is obsolete.
   */
   bool MdPairPotential::isPairListCurrent() 
   { return pairList_.isCurrent(boundary()); }

}
#endif
#endif
