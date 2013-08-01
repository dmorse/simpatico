#ifndef MCMD_STRUCTURE_FACTOR_CPP
#define MCMD_STRUCTURE_FACTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "StructureFactor.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/McMd_mpi.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/math/Constants.h>
#include <util/space/Dimension.h>
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   StructureFactor::StructureFactor(System& system) 
    : SystemDiagnostic<System>(system),
      isFirstStep_(true),
      isInitialized_(false)
   {  setClassName("StructureFactor"); }

   /*
   * Destructor.
   */
   StructureFactor::~StructureFactor() 
   {}

   /*
   * Read parameters from file, and allocate memory.
   */
   void StructureFactor::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "nMode", nMode_);
      nAtomType_ = system().simulation().nAtomType();
      modes_.allocate(nMode_, nAtomType_);
      readDMatrix<double>(in, "modes", modes_, nMode_, nAtomType_);
      read<int>(in, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, nWave_);

      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      structureFactors_.allocate(nWave_, nMode_);

      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void StructureFactor::loadParameters(Serializable::IArchive& ar)
   {
      Diagnostic::loadParameters(ar);
      ar & nAtomType_;
      loadParameter<int>(ar, "nMode", nMode_);
      loadDMatrix<double>(ar, "modes", modes_, nMode_, nAtomType_);
      loadParameter<int>(ar, "nWave", nWave_);
      loadDArray<IntVector>(ar, "waveIntVectors", waveIntVectors_, nWave_);
      ar & structureFactors_;
      ar & nSample_;

      // Validate
      if (nAtomType_ != system().simulation().nAtomType()) {
         UTIL_THROW("Inconsistent values of nAtomType_");
      }
      if (modes_.capacity1() != nMode_) {
         UTIL_THROW("Inconsistent capacity1 for modes array");
      }
      if (modes_.capacity2() != nAtomType_) {
         UTIL_THROW("Inconsistent capacity2 for modes array");
      }
      if (waveIntVectors_.capacity() != nWave_) {
         UTIL_THROW("Inconsistent capacity for waveIntVector");
      }

      // Allocate temporary data structures
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void StructureFactor::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /*
   * Clear accumulators.
   */
   void StructureFactor::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
      assert (nWave_ > 0);
      assert (nMode_ > 0);

      // Clear accumulators
      int i, j;
      for (i = 0; i < nWave_; ++i) {
         for (j = 0; j < nMode_; ++j) {
            structureFactors_(i, j) = 0.0;
         }
      }
      nSample_ = 0;
   }

   /* 
   * Increment structure factors for all wavevectors and modes.
   */
   void StructureFactor::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         fileMaster().openOutputFile(outputFileName("_max.dat"), outputFile_, !isFirstStep_);
         isFirstStep_ = false;

         Vector  position;
         std::complex<double>  expFactor;
         double  product;
         System::ConstMoleculeIterator  molIter;
         Molecule::ConstAtomIterator  atomIter;
         int  nSpecies, iSpecies, typeId, i, j;

         makeWaveVectors();

         // Set all Fourier modes to zero
         for (i = 0; i < nWave_; ++i) {
            for (j = 0; j < nMode_; ++j) {
               fourierModes_(i, j) = std::complex<double>(0.0, 0.0);
            }
         }
 
         // Loop over all atoms
         nSpecies = system().simulation().nSpecies();
         for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            system().begin(iSpecies, molIter); 
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(atomIter); 
               for ( ; atomIter.notEnd(); ++atomIter) {
                  position = atomIter->position();
                  typeId   = atomIter->typeId();
 
		  // Loop over wavevectors
                  for (i = 0; i < nWave_; ++i) {
               
                     product = position.dot(waveVectors_[i]);
                     expFactor = exp( product*Constants::Im );
                     for (j = 0; j < nMode_; ++j) {
                        fourierModes_(i, j) += modes_(j, typeId)*expFactor;
                     }
		 
                  }
  
               }
            }

         }

         // Increment structure factors
         double volume = system().boundary().volume();
         double norm;
         for (j = 0; j < nMode_; ++j) {
            double maxValue = 0.0;
            IntVector maxIntVector;
            double maxQ;
            for (i = 0; i < nWave_; ++i) {
               norm = std::norm(fourierModes_(i, j));
               if (double(norm/volume) >= maxValue) {
                  maxValue = double(norm/volume);
                  maxIntVector = waveIntVectors_[i];
                  maxQ = waveVectors_[i].abs();
               }
               structureFactors_(i, j) += norm/volume;
            }

            // Output current maximum S(q)
            outputFile_ << maxIntVector;
            outputFile_ << Dbl(maxQ,20,8);
            outputFile_ << Dbl(maxValue,20,8);
            outputFile_ << std::endl;
         }

         ++nSample_;

         outputFile_ << std::endl;
         outputFile_.close();
      }

   }

   /*
   * Calculate floating point wavevectors, using current boundary.
   */
   void StructureFactor::makeWaveVectors() 
   {
      Boundary* boundaryPtr = &system().boundary();
      Vector  dWave;
      int  i, j;
      for (i = 0; i < nWave_; ++i) {
         waveVectors_[i] = Vector::Zero;
         for (j = 0; j < Dimension; ++j) {
            dWave  = boundaryPtr->reciprocalBasisVector(j);
            dWave *= waveIntVectors_[i][j];
            waveVectors_[i] += dWave;
         }
      }
   }

   /*
   * Output final results to output file.
   */
   void StructureFactor::output() 
   {
      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      // Output structure factors to one file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      double      value;
      int         i, j, k;
      for (i = 0; i < nWave_; ++i) {
         for (k = 0; k < Dimension; ++k) {
            outputFile_ << Int(waveIntVectors_[i][k], 5);
         }
         outputFile_ << Dbl(waveVectors_[i].abs(), 20, 8);
         for (j = 0; j < nMode_; ++j) {
            value = structureFactors_(i, j)/double(nSample_);
            outputFile_ << Dbl(value, 18, 8);
         }
         outputFile_ << std::endl;
      }
      outputFile_.close();

      #if 0
      // Output each structure factor to a separate file
      std::string suffix;
      int         typeId;
      for (j = 0; j < nMode_; ++j) {

         // Construct file suffix for this structure factor
         suffix = std::string(".");
         suffix += toString(j);
         suffix += std::string(".dat");

         fileMaster().openOutputFile(outputFileName(suffix), outputFile_);

         // Loop over waves to output structure factor
         for (i = 0; i < nWave_; ++i) {
            for (k = 0; k < Dimension; ++k) {
               outputFile_ << Int(waveIntVectors_[i][k], 5);
            }
            outputFile_ << Dbl(waveVectors_[i].abs(), 20, 8);
            value = structureFactors_(i, j)/double(nSample_);
            outputFile_ << Dbl(value, 20, 8);
            outputFile_ << std::endl;
         }

         outputFile_.close();
      }
      #endif

   }

}
#endif
