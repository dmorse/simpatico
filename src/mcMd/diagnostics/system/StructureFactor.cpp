#ifndef MCMD_STRUCTURE_FACTOR_CPP
#define MCMD_STRUCTURE_FACTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
#include <mcMd/util/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/util/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   StructureFactor::StructureFactor(System& system) 
    : SystemDiagnostic<System>(system),
      isInitialized_(false)
   {}

   StructureFactor::~StructureFactor() 
   {}

   /// Read parameters from file, and allocate data array.
   void StructureFactor::readParam(std::istream& in) 
   {

      // Read interval and parameters for AutoCorrArray
      //SystemDiagnostic<System>::readParam(in);
      readInterval(in);
      readOutputFileName(in);

      nAtomType_ = system().simulation().nAtomType();

      read<int>(in, "nMode", nMode_);
      modes_.allocate(nMode_, nAtomType_);
      readDMatrix<double>(in, "modes", modes_, nMode_, nAtomType_);

      read<int>(in, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, nWave_);

      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      structureFactors_.allocate(nWave_, nMode_);

      maximumValue_.allocate(Samples);
      maximumWaveIntVector_.allocate(Samples);
      maximumQ_.allocate(Samples);

      isInitialized_ = true;
   }

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

      for (i=0; i < Samples; ++i) {
         maximumValue_[i] = 0.0;
      }

      nSample_ = 0;
   }
 
   /// Increment Structure Factor
   void StructureFactor::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

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
         for (i = 0; i < nWave_; ++i) {
            for (j = 0; j < nMode_; ++j) {
               norm = std::norm(fourierModes_(i, j));
               structureFactors_(i, j) += norm/volume;
               if (structureFactors_(i,j) >= maximumValue_[nSample_]) {
                  maximumValue_[nSample_] = structureFactors_(i,j);
                  maximumWaveIntVector_[nSample_] = waveIntVectors_[i];
                  maximumQ_[nSample_] = waveVectors_[i].abs();
               }
            }
         }

         ++nSample_;

      }

   }

   /**
   * Calculate floating point wavevectors.
   */
   void StructureFactor::makeWaveVectors() 
   {
      Vector    dWave;
      Boundary* boundaryPtr = &system().boundary();
      int       i, j;

      // Calculate wavevectors
      for (i = 0; i < nWave_; ++i) {
         waveVectors_[i] = Vector::Zero;
         for (j = 0; j < Dimension; ++j) {
            dWave  = boundaryPtr->reciprocalBasisVector(j);
            dWave *= waveIntVectors_[i][j];
            waveVectors_[i] += dWave;
            //std::cout << waveVectors_[i] << std::endl;
         }
      }
   }

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

      // Outputs history of maximum structure factors
      fileMaster().openOutputFile(outputFileName("_max.dat"), outputFile_);
      for (int i = 0; i < nSample_; ++i) {
         outputFile_ << maximumWaveIntVector_[i];
         outputFile_ << Dbl(maximumQ_[i], 20, 8);
         outputFile_ << Dbl(maximumValue_[i], 20, 8);
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

   /*
   * Save state to binary file archive.
   */
   void StructureFactor::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void StructureFactor::load(Serializable::IArchiveType& ar)
   { ar & *this; }
}
#endif
