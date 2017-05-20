#ifndef SIMP_NOPAIR
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McPairEnergyAverage.h"        // class header

#include <util/misc/FileMaster.h>  
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <util/archives/Serializable_includes.h>


#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   McPairEnergyAverage::McPairEnergyAverage(McSystem& system)
    : SystemAnalyzer<McSystem>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {
      setClassName("McPairEnergyAverage"); 
      selector_.setAvoidDoubleCounting(true); 
   }


   /*
   * Read parameters and initialize.
   */
   void McPairEnergyAverage::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      read<PairSelector>(in,"selector", selector_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (nSamplePerBlock_) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void McPairEnergyAverage::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);  
      loadParameter<int>(ar,"nSamplePerBlock", nSamplePerBlock_);
      loadParameter<PairSelector>(ar,"selector", selector_);
      ar & accumulator_;
 
      if (nSamplePerBlock_ != accumulator_.nSamplePerBlock()) {
         UTIL_THROW("Inconsistent values of nSamplePerBlock");
      } 

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (nSamplePerBlock_) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void McPairEnergyAverage::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /*
   * Clear accumulator.
   */
   void McPairEnergyAverage::setup() 
   {  accumulator_.clear(); }
 
   /*
   * Evaluate energy per particle, and add to ensemble. 
   */
   void McPairEnergyAverage::sample(long iStep) 
   {

      if (!isAtInterval(iStep)) return;

      Atom  *iAtomPtr, *jAtomPtr;
      int    nNeighbor, nInCell;
      int    ic, nc, i, j;
      double energy;
      double rsq;

      // Loop over cells
      energy = 0.0;
      nc = system().pairPotential().cellList().totCells();
      for (ic = 0; ic < nc; ++ic) {

         // Get array of neighbors_
         system().pairPotential().cellList().getCellNeighbors(ic, neighbors_, nInCell);
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

                     rsq = system().boundary().
                           distanceSq(iAtomPtr->position(), jAtomPtr->position());
                     energy += system().pairPotential().
                               energy(rsq, iAtomPtr->typeId(), jAtomPtr->typeId());

                  }

               }

            } // secondary atoms

         } // primary atoms

      } // cells

      accumulator_.sample(energy, outputFile_);
   }

   /*
   * Output results to file after simulation is completed.
   */
   void McPairEnergyAverage::output() 
   { 
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }

      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }
   
}
#endif 
