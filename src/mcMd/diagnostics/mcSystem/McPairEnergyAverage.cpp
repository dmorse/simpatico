#ifndef INTER_NOPAIR
#ifndef MCMD_MC_PAIR_ENERGY_AVERAGE_CPP
#define MCMD_MC_PAIR_ENERGY_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McPairEnergyAverage.h"        // class header

#include <mcMd/util/FileMaster.h>  
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
    : SystemDiagnostic<McSystem>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {
      selector_.setAvoidDoubleCounting(true); 
   }


   /*
   * Read parameters and initialize.
   */
   void McPairEnergyAverage::readParam(std::istream& in)
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
      int    ic, nc, i, j, iId, jId, iType, jType;
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
            iId      = iAtomPtr->id();
            iType    = iAtomPtr->typeId();
          
            // Loop over secondary atoms in this and neighboring cells
            for (j = 0; j < nNeighbor; ++j) {
               jAtomPtr = neighbors_[j];
               jId      = jAtomPtr->id();
               jType    = jAtomPtr->typeId();
     
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
   
   /*
   * Save state to binary file archive.
   */
   void McPairEnergyAverage::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void McPairEnergyAverage::load(Serializable::IArchiveType& ar)
   { ar & *this; }

}
#endif
#endif 
