/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomMSD.h"
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/archives/Serializable_includes.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   AtomMSD::AtomMSD(System& system) 
    : SystemAnalyzer<System>(system),
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
      if (speciesId_ >= system().simulation().nSpecies()) 
         UTIL_THROW("speciesId > nSpecies");
      if (atomId_ < 0)
         UTIL_THROW("Negative atomId");
      if (capacity_ <= 0)       
         UTIL_THROW("Negative capacity");

      Species* speciesPtr = &system().simulation().species(speciesId_);

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
   * Load state from an archive.
   */
   void AtomMSD::loadParameters(Serializable::IArchive& ar)
   { 
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "atomId", atomId_);
      loadParameter<int>(ar, "capacity", capacity_);

      // Validate parameters
      if (speciesId_ < 0)
         UTIL_THROW("Negative speciesId");
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      if (atomId_ < 0) UTIL_THROW("Negative atomId");
      if (capacity_ <= 0) UTIL_THROW("Negative capacity");

      Species* speciesPtr = &system().simulation().species(speciesId_);

      if (atomId_ >= speciesPtr->nAtom()) {
         UTIL_THROW("atomId >= nAtom");
      }

      // Maximum possible number of molecules of this species
      int speciesCapacity = speciesPtr->capacity();

      // Allocate local arrays
      truePositions_.allocate(speciesCapacity);
      oldPositions_.allocate(speciesCapacity);
      shifts_.allocate(speciesCapacity);

      // Allocate memory for the AutoCorrArray accumulator
      accumulator_.setParam(speciesCapacity, capacity_);

      ar & accumulator_;
      ar & truePositions_;
      ar & oldPositions_;
      ar & shifts_;
      ar & nMolecule_;

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void AtomMSD::save(Serializable::OArchive& ar)
   { ar & (*this); }

   /*
   * Initialize at beginning of simulation.
   */
   void AtomMSD::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }

      // Set number of molecules of this species in the System. 
      nMolecule_ = system().nMolecule(speciesId_);
      accumulator_.setNEnsemble(nMolecule_);
      accumulator_.clear();

      // Store initial positions, and set initial shift vectors.
      Vector r;
      IntVector  zero(0);
      for (int i = 0; i < nMolecule_; ++i) {
         r = system().molecule(speciesId_, i).atom(atomId_).position();
         system().boundary().shift(r);
         oldPositions_[i] = r;
         shifts_[i] = zero;
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
      if (nMolecule_ != system().nMolecule(speciesId_)) {
         UTIL_THROW("Number of molecules has changed.");
      }

      Vector     r;
      IntVector  shift;
      Vector     lengths = system().boundary().lengths();
      int        i, j;
      for (i = 0; i < nMolecule_; ++i) {

         r = system().molecule(speciesId_, i).atom(atomId_).position();
         system().boundary().shift(r);

         // Compare current r to previous position, oldPositions_[i]
         system().boundary().distanceSq(r, oldPositions_[i], shift);

         // If this atom crossed a boundary, increment its shift vector
         shifts_[i] += shift;

         // Reconstruct true position
         for (j = 0; j < Dimension; ++j) {
            truePositions_[i][j] = r[j] + shifts_[i][j]*lengths[j];
         }

         // Store current position in box for comparison to next one
         oldPositions_[i] = r;
  
      }
      accumulator_.sample(truePositions_);

   }

   /// Output results to file after simulation is completed.
   void AtomMSD::output() 
   {  

      // Output parameters
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
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
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}
