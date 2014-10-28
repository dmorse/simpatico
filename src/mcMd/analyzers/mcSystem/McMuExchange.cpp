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

#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
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

         #if 0
         Species* speciesPtr;
         Molecule* molPtr;
         Atom* atom1Ptr;
         Atom* atom2Ptr;

         double beta = energyEnsemble().beta();
         speciesPtr = &(simulation().species(speciesId_));

         // Loop over molecules in species {
         //   for (int j=0; flipAtomIds__.size(); ++j) {
         //      i = flipAtomIds__[j];
         //      Get neighbors from cell list.
         //      Calculate old and new pair energies.
         //   }
         // } 
         #endif

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
