#ifndef SPAN_ATOM_MSD_CPP
#define SPAN_ATOM_MSD_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomMSD.h"
#include <spAn/chemistry/Molecule.h>
#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Species.h>
#include <spAn/processor/Processor.h>
#include <util/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/archives/Serializable_includes.h>

#include <util/global.h>

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   AtomMSD::AtomMSD(Processor& configuration) 
    : Analyzer(configuration),
      outputFile_(),
      accumulator_(),
      truePositions_(),
      oldPositions_(),
      shifts_(),
      speciesId_(-1),
      atomId_(-1),
      nMolecule_(-1),
      capacity_(-1),
      isInitialized_(false)
   {  setClassName("AtomMSD"); }

   /*
   * Read parameters from file, and allocate data array.
   */
   void AtomMSD::readParameters(std::istream& in) 
   {

      // Read interval and parameters for AutoCorrArray
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "atomId", atomId_);
      read<int>(in, "capacity", capacity_);

      // Validate parameters
      if (speciesId_ < 0)
         UTIL_THROW("Negative speciesId");
      if (speciesId_ >= processor().nSpecies()) 
         UTIL_THROW("speciesId > nSpecies");
      if (atomId_ < 0)
         UTIL_THROW("Negative atomId");
      if (capacity_ <= 0)       
         UTIL_THROW("Negative capacity");

      Species* speciesPtr = &processor().species(speciesId_);

      if (atomId_ >= speciesPtr->nAtom()) 
         UTIL_THROW("atomId >= nAtom");

      // Maximum possible number of molecules of this species
      int speciesCapacity = speciesPtr->capacity();

      // Allocate local arrays
      truePositions_.allocate(speciesCapacity);
      oldPositions_.allocate(speciesCapacity);
      shifts_.allocate(speciesCapacity);

      // Allocate memory for the AutoCorrArray accumulator
      accumulator_.setParam(speciesCapacity, capacity_);

      isInitialized_ = true;
   }

   /*
   * Initialize at beginning of simulation.
   */
   void AtomMSD::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }

      // Set number of molecules of this species in the Processor. 
      nMolecule_ = processor().species(speciesId_).size();
      accumulator_.setNEnsemble(nMolecule_);
      accumulator_.clear();

      // Store initial positions, and set initial shift vectors.
      Vector r;
      IntVector zero(0);
      Species::MoleculeIterator iter;
      int i = 0;

      processor().species(speciesId_).begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         r = iter->atom(atomId_).position;
         processor().boundary().shift(r);
         oldPositions_[i] = r;
         shifts_[i] = zero;
         ++i;
      }
   }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void AtomMSD::sample(long iStep) 
   { 
      if (!isAtInterval(iStep)) return;

      // Confirm that nMolecule has remained constant, and nMolecule > 0.
      if (nMolecule_ <= 0) {
         UTIL_THROW("nMolecule <= 0");
      }
      if (nMolecule_ != processor().species(speciesId_).size()) {
         UTIL_THROW("Number of molecules has changed.");
      }

      Vector r;
      IntVector shift;
      Vector lengths = processor().boundary().lengths();
      int i, j;
      Species::MoleculeIterator iter;

      i = 0;
      processor().species(speciesId_).begin(iter);
      for ( ; iter.notEnd(); ++iter) {

         r = iter->atom(atomId_).position;
         processor().boundary().shift(r);

         // Compare current r to previous position, oldPositions_[i]
         processor().boundary().distanceSq(r, oldPositions_[i], shift);

         // If this atom crossed a boundary, increment its shift vector
         shifts_[i] += shift;

         // Reconstruct true position
         for (j = 0; j < Dimension; ++j) {
            truePositions_[i][j] = r[j] + shifts_[i][j]*lengths[j];
         }

         // Store current position in box for comparison to next one
         oldPositions_[i] = r;

         ++i;  
      }
      accumulator_.sample(truePositions_);

   }

   #if 1
   /// Output results to file after simulation is completed.
   void AtomMSD::output() 
   {  

      // Output parameters
      if (processor().hasFileMaster()) {
         processor().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      } else {
         outputFile_.open(outputFileName(".prm").c_str());
      }

      writeParam(outputFile_); 
      outputFile_ << std::endl;
      outputFile_ << std::endl;
      outputFile_ << "nMolecule      " << accumulator_.nEnsemble() 
                  << std::endl;
      outputFile_ << "buffercapacity " << accumulator_.bufferCapacity()  
                  << std::endl;
      outputFile_ << "nSample        " << accumulator_.nSample() 
                  << std::endl;
      outputFile_ << std::endl;
      //outputFile_ << "Format of *.dat file:" << std::endl;
      //outputFile_ << "[i]  [MSD]" << std::endl;
      outputFile_.close();

      // Output statistical analysis to separate data file
      if (processor().hasFileMaster()) {
         processor().fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      } else {
         outputFile_.open(outputFileName(".dat").c_str());
      }
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }
   #endif

}
#endif 
