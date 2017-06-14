/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
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
      SystemInterface(system)
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
      // Precondition
      if (!pairList_.isInitialized()) {
         UTIL_THROW("PairList not initialized in MdPairPotential::buildPairList");
      }

      // Set up an empty PairList with an empty internal CellList.
      pairList_.setup(boundary());

      // Add every Atom in this System to the CellList
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
