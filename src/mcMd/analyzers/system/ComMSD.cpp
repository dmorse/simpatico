/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ComMSD.h"
#include <mcMd/simulation/System.h>               // base class template parameter
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>                 // forward declaration in ComMSD.h
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
   ComMSD::ComMSD(System& system) 
    : SystemAnalyzer<System>(system),
      outputFile_(),
      accumulator_(),
      truePositions_(),
      oldPositions_(),
      shifts_(),
      speciesId_(-1),
      nMolecule_(-1),
      nAtom_(-1),
      capacity_(-1),
      isInitialized_(false)
   { setClassName("ComMSD"); }

   /*
   * Read parameters from file, and allocate data array.
   */
   void ComMSD::readParameters(std::istream& in) 
   {
      // Read interval and parameters for AutoCorrArray
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "capacity", capacity_);

      // Validate parameters
      if (speciesId_ < 0)       UTIL_THROW("Negative speciesId");
      if (speciesId_ >= system().simulation().nSpecies()) 
                                UTIL_THROW("speciesId > nSpecies");
      if (capacity_ <= 0)       UTIL_THROW("Negative capacity");

      // Maximum possible number of molecules of this species
      int speciesCapacity;
      speciesCapacity = system().simulation().species(speciesId_).capacity();

      // Allocate local arrays
      truePositions_.allocate(speciesCapacity);
      oldPositions_.allocate(speciesCapacity);
      shifts_.allocate(speciesCapacity);

      // Allocate memory for the AutoCorrArray accumulator
      accumulator_.setParam(speciesCapacity, capacity_);

      // Get number of atoms per molecule.
      nAtom_ = system().simulation().species(speciesId_).nAtom();

      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void ComMSD::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "capacity", capacity_);
      ar & nAtom_;

      // Validate parameters
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      if (capacity_ <= 0) {
         UTIL_THROW("Negative capacity");
      }
      if (nAtom_ != system().simulation().species(speciesId_).nAtom()) {
         UTIL_THROW("Inconsistent values for nAtom");
      }

      // Maximum possible number of molecules of this species
      int speciesCapacity;
      speciesCapacity = system().simulation().species(speciesId_).capacity();

      // Allocate memory
      truePositions_.allocate(speciesCapacity);
      oldPositions_.allocate(speciesCapacity);
      shifts_.allocate(speciesCapacity);
      accumulator_.setParam(speciesCapacity, capacity_);

      // Finish reading internal state
      ar & accumulator_;
      ar & truePositions_;
      ar & oldPositions_;
      ar & shifts_;
      ar & nMolecule_;
   }

   /*
   * Save state to archive.
   */
   void ComMSD::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Set number of molecules and initialize state.
   */
   void ComMSD::setup() 
   {
      // Precondition: Confirm that readParam was called
      if (!isInitialized_)
         UTIL_THROW("Analyzer not initialized");

      // Get number of molecules of this species in this System. 
      nMolecule_ = system().nMolecule(speciesId_);
      if (nMolecule_ <= 0)
         UTIL_THROW("nMolecule_ <= 0");
      accumulator_.setNEnsemble(nMolecule_);

      accumulator_.clear();

      // Store initial positions, and set initial shift vectors.
      Vector r;
      for (int i = 0; i < nMolecule_; ++i) {
         r = system().molecule(speciesId_, i).atom(0).position();
         system().boundary().shift(r);
         oldPositions_[i] = r;
         shifts_[i] = IntVector(0);
      }

   }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void ComMSD::sample(long iStep) 
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

         r = system().molecule(speciesId_, i).atom(0).position();
         system().boundary().shift(r);

         // Compare current r to previous position, oldPositions_[i]
         system().boundary().distanceSq(r, oldPositions_[i], shift);

         // If atom 0 crossed a boundary, increment its shift vector
         shifts_[i] += shift;

         // Store current position in box for comparison to next one
         oldPositions_[i] = r;
 
         // Add the other particle's positions to C.O.M.
         moleculeCom(i, r);

         // Reconstruct true position
         for (j = 0; j < Dimension; ++j) {
            truePositions_[i][j] = r[j] + double(shifts_[i][j])*lengths[j];
         }
 
      }
      accumulator_.sample(truePositions_);

   }

   /*
   * Output results to file after simulation is completed.
   */
   void ComMSD::output() 
   {  

      // Echo parameter to analyzer log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_ << std::endl;
      outputFile_ << std::endl;
      outputFile_ << "nMolecule      " << accumulator_.nEnsemble() << std::endl;
      outputFile_ << "bufferCapacity " << accumulator_.bufferCapacity() 
                  << std::endl;
      outputFile_ << "nSample        " << accumulator_.nSample()   << std::endl;
      outputFile_ << std::endl;
      outputFile_ << "Format of *.dat file:" << std::endl;
      outputFile_ << "[i]  [MSD]"
                  << std::endl;
      outputFile_.close();

      // Output statistical analysis to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();
   }

   /*
   * Evaluate center of mass of a single molecule, given a starting vector R0.
   */
   void ComMSD::moleculeCom(int iMol, Vector &R0)
   {
      Vector R = R0, r1, r2, dr;

      // Assuming a linear chain (including Ring molecules).
      for (int i = 1; i < nAtom_-1; ++i) {
         r1 = system().molecule(speciesId_, iMol).atom(i-1).position();
         r2 = system().molecule(speciesId_, iMol).atom(i).position();
         system().boundary().distanceSq(r2, r1, dr);
         R  += dr;
         R0 += R;
      }

      // Do the average.
      R0 /= double(nAtom_);
   }

}
