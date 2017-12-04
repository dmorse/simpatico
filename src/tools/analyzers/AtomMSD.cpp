/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomMSD.h"
#include <tools/chemistry/Molecule.h>
#include <tools/chemistry/Atom.h>
#include <tools/chemistry/Species.h>
#include <tools/storage/Configuration.h>
#include <simp/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/archives/Serializable_includes.h>

#include <util/global.h>

namespace Tools
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   AtomMSD::AtomMSD(Processor& processor) 
    : Analyzer(processor),
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
   * Constructor.
   */
   AtomMSD::AtomMSD(Configuration& configuration, FileMaster& fileMaster) 
    : Analyzer(configuration, fileMaster),
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
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= configuration().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      Species& species = configuration().species(speciesId_);
      if (atomId_ < 0) {
         UTIL_THROW("Negative atomId");
      }
      if (atomId_ >= species.nAtom()) {
         UTIL_THROW("atomId >= nAtom");
      }
      if (capacity_ <= 0) {
         UTIL_THROW("Negative capacity");
      }

      // Maximum possible number of molecules of this species
      int speciesCapacity = species.capacity();

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
      Species& species = configuration().species(speciesId_);
      Boundary& boundary = configuration().boundary();

      // Set number of molecules of this species in the Configuration. 
      nMolecule_ = species.size();
      accumulator_.setNEnsemble(nMolecule_);
      accumulator_.clear();

      // Store initial positions, and set initial shift vectors.
      Vector r;
      IntVector zero(0);
      Species::Iterator iter;
      int i = 0;

      species.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         r = iter->atom(atomId_).position;
         boundary.shift(r);
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
      if (isAtInterval(iStep)) {
         Species&  species  = configuration().species(speciesId_);
         Boundary& boundary = configuration().boundary();
   
         // Confirm that nMolecule is positive, and unchanged.
         if (nMolecule_ <= 0) {
            UTIL_THROW("nMolecule <= 0");
         }
         if (nMolecule_ != species.size()) {
            UTIL_THROW("Number of molecules has changed.");
         }
   
         Vector r;
         IntVector shift;
         Vector lengths = boundary.lengths();
         int i, j;
         Species::Iterator iter;
   
         i = 0;
         for (species.begin(iter); iter.notEnd(); ++iter) {
   
            r = iter->atom(atomId_).position;
            boundary.shift(r);
   
            // Compare current r to previous position, oldPositions_[i]
            boundary.distanceSq(r, oldPositions_[i], shift);
   
            // If a boundary was crossed, increment shift vector
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

   }

   /// Output results to file after simulation is completed.
   void AtomMSD::output() 
   {  

      // Output parameters
      fileMaster().openOutputFile(outputFileName(".prm"), 
                                              outputFile_);

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
      fileMaster().openOutputFile(outputFileName(".dat"), 
                                              outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}
