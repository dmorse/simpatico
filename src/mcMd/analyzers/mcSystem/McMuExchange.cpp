#ifndef MCMD_MC_MU_EXCHANGE_CPP
#define MCMD_MC_MU_EXCHANGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMuExchange.h"                            // class header
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>

#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/neighbor/CellList.h>

#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <util/space/Vector.h>
#include <util/global.h>

#include <math.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McMuExchange::McMuExchange(McSystem& system)
    : SystemAnalyzer<McSystem>(system),
      simulationPtr_(&system.simulation()),
      boundaryPtr_(&system.boundary()),
      outputFile_(),
      accumulator_(),
      newTypeIds_(),
      oldTypeIds_(),
      flipAtomIds_(),
      speciesId_(-1),
      nAtom_(-1),
      isInitialized_(false)
   {}


   /*
   * Read parameters and initialize.
   */
   void McMuExchange::readParameters(std::istream& in)
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId >= nSpecies");
      }
      Species* speciesPtr = &(system().simulation().species(speciesId_));
      nAtom_ = speciesPtr->nAtom();
      newTypeIds_.allocate(nAtom_);
      readDArray<int>(in, "newTypeIds", newTypeIds_, nAtom_);

      oldTypeIds_.allocate(nAtom_);
      flipAtomIds_.allocate(nAtom_);
      for (int i = 0; i < nAtom_; ++i) {
         oldTypeIds_[i] = speciesPtr->atomTypeId(i);
         if (newTypeIds_[i] != oldTypeIds_[i]) {
            flipAtomIds_.append(i);
         }
      }
      flipAtomIds_.allocate(nAtom_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

      isInitialized_ = true;
   }

   /*
   * Clear accumulator.
   */
   void McMuExchange::setup()
   {  accumulator_.clear(); }

   /*
   * Evaluate change in energy, add Boltzmann factor to accumulator.
   */
   void McMuExchange::sample(long iStep)
   {
      if (isAtInterval(iStep))  {
         if (!system().energyEnsemble().isIsothermal()) {
            UTIL_THROW("EnergyEnsemble is not isothermal");
         }

         //Species& species = simulation().species(speciesId_);
         McPairPotential& potential = system().pairPotential();
         const CellList& cellList = potential.cellList();
         Atom* ptr1 = 0;
         Atom* ptr2 = 0;
         Mask* maskPtr = 0;
         double beta = system().energyEnsemble().beta();
         double rsq, dE, boltzmann;
         int i, j, k, nNeighbor;
         int id1, id2, t1, t2, t1New;


         System::MoleculeIterator molIter;
         Molecule::AtomIterator atomIter;
         system().begin(speciesId_, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            dE = 0.0;
            for (j=0; flipAtomIds_.size(); ++j) {
               i = flipAtomIds_[j];
               t1New = newTypeIds_[i];
               ptr1 = &molIter->atom(i);
               t1 = ptr1->typeId();
               id1 = ptr1->id();
               maskPtr = &(ptr1->mask());
               cellList.getNeighbors(ptr1->position(), neighbors_);
               nNeighbor = neighbors_.size();
               for (k = 0; k < nNeighbor; ++k) {
                  ptr2 = neighbors_[k];
                  t2 = ptr2->typeId();
                  id2 = ptr2->id();

                  // Check if atoms are identical
                  if (id1 != id2) {
          
                     // Check if pair is masked
                     if (!maskPtr->isMasked(*ptr2)) {
                        rsq = boundary().distanceSq(ptr1->position(), 
                                                    ptr2->position());
                        dE += potential.energy(rsq, t1New, t2);
                        dE -= potential.energy(rsq, t1, t2);
                     }
                  }
               }
            }
            boltzmann = exp(-beta*dE);
            accumulator_.sample(boltzmann);
         }
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void McMuExchange::output()
   {
      #if 0
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }
      #endif

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
   void McMuExchange::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void McMuExchange::load(Serializable::IArchive& ar)
   { ar & *this; }

}
#endif
