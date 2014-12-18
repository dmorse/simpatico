/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PairEnergy.h"
#include <tools/chemistry/Atom.h>
#include <tools/storage/Configuration.h>
#include <util/space/Vector.h>
#include <tools/neighbor/CellList.h>
#include <tools/neighbor/Cell.h>

namespace Tools
{

   /*
   * Constructor.
   */
   PairEnergy::PairEnergy(Configuration& configuration)
    : Analyzer(configuration),
      interaction_(),
      cellList_(),
      cutoff_(0.0),
      atomCapacity_(0),
      isInitialized_(false)
   {  setClassName("PairEnergy"); }

   /*
   * Constructor.
   */
   PairEnergy::PairEnergy(Processor& processor)
    : Analyzer(processor),
      interaction_(),
      cellList_(),
      cutoff_(0.0),
      atomCapacity_(0),
      isInitialized_(false)
   {  setClassName("PairEnergy"); }

   /*
   * Constructor.
   */
   PairEnergy::PairEnergy(Configuration& configuration, FileMaster& fileMaster)
    : Analyzer(configuration, fileMaster),
      interaction_(),
      cellList_(),
      cutoff_(0.0),
      atomCapacity_(0),
      isInitialized_(false)
   {  setClassName("PairEnergy"); }


   void PairEnergy::readParameters(std::istream& in) 
   {
      // Read interval and parameters for AutoCorrArray
      readInterval(in);
      readOutputFileName(in);

      read<int>(in, "atomCapacity", atomCapacity_);
      interaction_.setNAtomType(1);
      interaction_.readParameters(in);
      
      cutoff_ = interaction_.maxPairCutoff();

      isInitialized_ = true;
   }


   /*
   * Initialize at beginning of simulation.
   */
   void PairEnergy::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }

      // setup cell List
      Vector lengths = configuration().boundary().lengths();
      
      Vector lower(0.0, 0.0, 0.0);
      Vector upper = lengths; 
      Vector cutoffs(cutoff_, cutoff_, cutoff_); 

      cellList_.allocate(atomCapacity_, lower, upper, cutoffs);
      cellList_.makeGrid(lower, upper, cutoffs);
   }

   void PairEnergy::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) return;

      // clear cell list
      cellList_.clear();

      AtomStorage::Iterator atomIter;
      for (configuration().atoms().begin(atomIter); atomIter.notEnd(); ++atomIter) {
         cellList_.placeAtom(*atomIter);
      }

      // build/update cell list
      cellList_.build();
      cellList_.update();

      if (!cellList_.isValid()) {
         UTIL_THROW("Cell List Invalid\n");
      }

      // Find all neighbor pairs within a cutoff (using cell list)
      Cell::NeighborArray neighbors;
      CellAtom* cellAtomPtr1;
      CellAtom* cellAtomPtr2;
      Vector dr;
      int   na = 0; // total number of atoms
      int   nn = 0; // number of neighbors in a cell
      int   np = 0; // Number of pairs within cutoff

      const Cell* cellPtr = cellList_.begin();
      double energy = 0.0;

      while (cellPtr) {
         cellPtr->getNeighbors(neighbors);
         na = cellPtr->nAtom();
         nn = neighbors.size();
         for (int i = 0; i < na; ++i) {
            cellAtomPtr1 = neighbors[i];
            for (int j = 0; j < na; ++j) {
               cellAtomPtr2 = neighbors[j];
               if (cellAtomPtr2 > cellAtomPtr1) {
                  double rsq = configuration().boundary().distanceSq(cellAtomPtr1->ptr()->position, cellAtomPtr2->ptr()->position, dr);
                  if (rsq <= cutoff_*cutoff_) {
                     ++np;
                     energy += interaction_.energy( rsq, cellAtomPtr1->ptr()->typeId, cellAtomPtr2->ptr()->typeId );
                  }
               }
            }
            for (int j = na; j < nn; ++j) {
               cellAtomPtr2 = neighbors[j];
               double rsq = configuration().boundary().distanceSq(cellAtomPtr1->ptr()->position, cellAtomPtr2->ptr()->position, dr);
               if (rsq <= cutoff_*cutoff_) {
                  ++np;
                  energy += interaction_.energy( rsq, cellAtomPtr1->ptr()->typeId, cellAtomPtr2->ptr()->typeId );
               }
            }
         }
         cellPtr = cellPtr->nextCellPtr();
      }

      timesteps_.append(iStep);
      energies_.append(energy);
   }

   /// Output results to file after simulation is completed.
   void PairEnergy::output() 
   {
      // Output parameters
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(), 
                                              outputFile_);
      char str[50];
      sprintf(&str[0], "%10s       %-7s", "TIMESTEP", "ENERGY");

      outputFile_ << str << std::endl;
      for (int i = 0; i < timesteps_.size(); i++) {
         sprintf(&str[0], "%10i       %-12.7e", timesteps_[i], energies_[i]);
         outputFile_ << str << std::endl;
      }

      outputFile_.close();
   }

}

