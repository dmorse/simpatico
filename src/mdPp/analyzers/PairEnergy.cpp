/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PairEnergy.h"
#include <mdPp/storage/Configuration.h>
#include <mdPp/neighbor/CellList.h>
#include <mdPp/neighbor/Cell.h>
#include <mdPp/chemistry/Atom.h>
#include <util/space/Vector.h>

namespace MdPp
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
   PairEnergy::PairEnergy(Configuration& configuration, 
                          FileMaster& fileMaster)
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

      // Setup cell List
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

      // Clear cell list
      cellList_.clear();

      // Place all atoms
      AtomStorage::Iterator atomIter;
      configuration().atoms().begin(atomIter); 
      for ( ; atomIter.notEnd(); ++atomIter) {
         cellList_.placeAtom(*atomIter);
      }

      // Build/update cell list
      cellList_.build();
      cellList_.update();

      if (!cellList_.isValid()) {
         UTIL_THROW("Cell List Invalid\n");
      }

      // Compute pair energy using loop over cells
      Vector dr;
      double energy = 0.0;
      double rsq;
      Cell::NeighborArray neighbors;
      Boundary& boundary = configuration().boundary();
      CellAtom* cellAtomPtr1 = 0;
      CellAtom* cellAtomPtr2 = 0;
      const Cell* cellPtr = 0;
      int na = 0; // total number of atoms
      int nn = 0; // number of neighbors in a cell
      int np = 0; // Number of pairs within cutoff

      // Loop over cells in CellList
      cellPtr = cellList_.begin();
      while (cellPtr) {
         cellPtr->getNeighbors(neighbors);
         na = cellPtr->nAtom();
         nn = neighbors.size();
         // Loop over primary atoms, from this cell
         for (int i = 0; i < na; ++i) {
            cellAtomPtr1 = neighbors[i];
            // Loop over secondary atoms in this cell
            for (int j = 0; j < na; ++j) {
               cellAtomPtr2 = neighbors[j];
               if (cellAtomPtr2 > cellAtomPtr1) {
                  rsq = boundary.distanceSq(cellAtomPtr1->ptr()->position, 
                                            cellAtomPtr2->ptr()->position, dr);
                  if (rsq <= cutoff_*cutoff_) {
                     ++np;
                     energy += interaction_.energy(rsq, 
                                            cellAtomPtr1->ptr()->typeId,
                                            cellAtomPtr2->ptr()->typeId);
                  }
               }
            }
            // Loop over atoms in neighboring cells
            for (int j = na; j < nn; ++j) {
               cellAtomPtr2 = neighbors[j];
               rsq = boundary.distanceSq(cellAtomPtr1->ptr()->position, 
                                         cellAtomPtr2->ptr()->position, dr);
               if (rsq <= cutoff_*cutoff_) {
                  ++np;
                  energy += interaction_.energy(rsq, 
                                         cellAtomPtr1->ptr()->typeId, 
                                         cellAtomPtr2->ptr()->typeId);
               }
            }
         }
         cellPtr = cellPtr->nextCellPtr();
      }

      timesteps_.append(iStep);
      energies_.append(energy);
   }

   /*
   * Output results to file after simulation is completed.
   */
   void PairEnergy::output() 
   {
      // Output parameters
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output energies
      fileMaster().openOutputFile(outputFileName(), outputFile_);
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

