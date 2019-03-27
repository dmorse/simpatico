#ifndef SIMP_NOPAIR
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McPairEnergyAnalyzer.h"        // class header

#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <simp/boundary/Boundary.h>
//#include <util/archives/Serializable_includes.h>


#include <cstdio> 

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /* 
   * Constructor.
   */
   McPairEnergyAnalyzer::McPairEnergyAnalyzer(McSystem& system)
    : AverageAnalyzer<McSystem>(system)
   {
      setClassName("McPairEnergyAnalyzer"); 
      selector_.setAvoidDoubleCounting(true); 
   }


   /*
   * Read parameters and initialize.
   */
   void McPairEnergyAnalyzer::readParameters(std::istream& in)
   {
      AverageAnalyzer<McSystem>::readParameters(in);
      read<PairSelector>(in, "selector", selector_);
   }

   /*
   * Load state from an archive.
   */
   void McPairEnergyAnalyzer::loadParameters(Serializable::IArchive& ar)
   {
      AverageAnalyzer<McSystem>::loadParameters(ar);  
      loadParameter<PairSelector>(ar, "selector", selector_);
   }

   /*
   * Save state to archive.
   */
   void McPairEnergyAnalyzer::save(Serializable::OArchive& ar)
   {
      AverageAnalyzer<McSystem>::save(ar);
      ar << selector_; 
   }

   /*
   * Evaluate energy per particle, and add to ensemble. 
   */
   void McPairEnergyAnalyzer::compute() 
   {
      McPairPotential& potential = system().pairPotential();
      Boundary& boundary = system().boundary();
      Atom  *iAtomPtr, *jAtomPtr;
      int nNeighbor, nInCell;
      int ic, nc, i, j;
      double rsq;

      // Loop over cells
      value_ = 0.0;
      nc = potential.cellList().totCells();
      for (ic = 0; ic < nc; ++ic) {

         // Get array of neighbors_
         potential.cellList().getCellNeighbors(ic, neighbors_, nInCell);
         nNeighbor = neighbors_.size();
  
         // Loop over primary atoms in this cell
         for (i = 0; i < nInCell; ++i) {
            iAtomPtr = neighbors_[i];
          
            // Loop over secondary atoms in this and neighboring cells
            for (j = 0; j < nNeighbor; ++j) {
               jAtomPtr = neighbors_[j];
     
               if (selector_.match(*iAtomPtr, *jAtomPtr)) {

                  // Exclude masked pairs
                  if (!iAtomPtr->mask().isMasked(*jAtomPtr)) {

                     rsq = boundary.distanceSq(iAtomPtr->position(), 
                                               jAtomPtr->position());
                     value_ += potential.energy(rsq, 
                               iAtomPtr->typeId(), jAtomPtr->typeId());

                  }

               }

            } // secondary atoms

         } // primary atoms

      } // cells
   }
   
}
#endif 
