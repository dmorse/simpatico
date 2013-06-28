#ifndef DDMD_STRUCTURE_FACTOR_CPP
#define DDMD_STRUCTURE_FACTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "StructureFactor.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/simulation/SimulationAccess.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <util/boundary/Boundary.h>
#include <util/math/Constants.h>
#include <util/space/Dimension.h>
#include <util/space/IntVector.h>
#include <util/misc/ioUtil.h>
#include <util/mpi/MpiLoader.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace DdMd
{

   using namespace Util;

   /// Constructor.
   StructureFactor::StructureFactor(Simulation& simulation) 
    : Diagnostic(simulation),
      isInitialized_(false)
   {  setClassName("StructureFactor"); }

   StructureFactor::~StructureFactor() 
   {}

   /// Read parameters from file, and allocate data array.
   void StructureFactor::readParameters(std::istream& in) 
   {
      nAtomType_ = simulation().nAtomType();

      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "nMode", nMode_);
      modes_.allocate(nMode_, nAtomType_);
      readDMatrix<double>(in, "modes", modes_, nMode_, nAtomType_);
      read<int>(in, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, nWave_);

      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      totalFourierModes_.allocate(nWave_, nMode_);

      if (simulation().domain().isMaster()) {

         structureFactors_.allocate(nWave_, nMode_);
         int i, j;
         for (i = 0; i < nWave_; ++i) {
            for (j = 0; j < nMode_; ++j) {
               structureFactors_(i, j) = 0.0;
            }
         }

         maximumValue_.allocate(nMode_);
         maximumWaveIntVector_.allocate(nMode_);
         maximumQ_.allocate(nMode_);
         for (j = 0; j < nMode_; ++j) {
            maximumValue_[j].reserve(Samples);
            maximumWaveIntVector_[j].reserve(Samples);
            maximumQ_[j].reserve(Samples);
         }

      }
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void StructureFactor::loadParameters(Serializable::IArchive &ar)
   {
      nAtomType_ = simulation().nAtomType();

      // Load and broadcast parameter file parameters
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar, "nMode", nMode_);
      modes_.allocate(nMode_, nAtomType_);
      loadDMatrix<double>(ar, "modes", modes_, nMode_, nAtomType_);
      loadParameter<int>(ar, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      loadDArray<IntVector>(ar, "waveIntVectors", waveIntVectors_, nWave_);

      // Load and broadcast nSample_
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);

      // Allocate and load accumulators that exist only on master.
      if (simulation().domain().isMaster()) {
         structureFactors_.allocate(nWave_, nMode_);
         ar >> structureFactors_;

         maximumValue_.allocate(nMode_);
         maximumWaveIntVector_.allocate(nMode_);
         maximumQ_.allocate(nMode_);
         for (int j = 0; j < nMode_; ++j) {
            ar >> maximumValue_[j];
            ar >> maximumWaveIntVector_[j];
            ar >> maximumQ_[j];
         }
      }

      // Allocate work space (all processors).
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      totalFourierModes_.allocate(nWave_, nMode_);

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void StructureFactor::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nMode_;
      ar << modes_;
      ar << nWave_;
      ar << waveIntVectors_;

      ar << nSample_;

      ar << structureFactors_;
      for (int j = 0; j < nMode_; ++j) {
         ar << maximumValue_[j];
         ar << maximumWaveIntVector_[j];
         ar << maximumQ_[j];
      }
   }
  
   /*
   * Clear accumulators.
   */
   void StructureFactor::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
      assert (nWave_ > 0);
      assert (nMode_ > 0);

      nSample_ = 0;
      if (simulation().domain().isMaster()) {

         int i, j;
         for (i = 0; i < nWave_; ++i) {
            for (j = 0; j < nMode_; ++j) {
               structureFactors_(i, j) = 0.0;
            }
         }
         for (int j = 0; j < nMode_; ++j) {
            maximumValue_[j].clear();
            maximumWaveIntVector_[j].clear();
            maximumQ_[j].clear();
         }
      }
   }

   /*
   * Increment structure factor.
   */
   void StructureFactor::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         Vector position;
         std::complex<double> expFactor;
         double  product;
         AtomIterator  atomIter;
         int i, j, typeId;

         makeWaveVectors();

         // Set all Fourier modes to zero
         for (i = 0; i < nWave_; ++i) {
            for (j = 0; j < nMode_; ++j) {
               fourierModes_(i, j) = std::complex<double>(0.0, 0.0);
            }
         }

         simulation().atomStorage().begin(atomIter);
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

         for (i = 0; i < nWave_; ++i) {
            for (j = 0; j < nMode_; ++j) {
               totalFourierModes_(i, j) = std::complex<double>(0.0, 0.0);
            }
         }
 
         #ifdef UTIL_MPI
         // Loop over wavevectors
         for (int i = 0; i < nWave_; ++i) {
            for (int j = 0; j < nMode_; ++j) {
            //Sum values from all processors.
            simulation().domain().communicator().
                         Reduce(&fourierModes_(i, j), &totalFourierModes_(i, j),
                                1, MPI::DOUBLE_COMPLEX, MPI::SUM, 0);
            }
         }
         #else
         for (int i = 0; i < nWave_; ++i) {
            for (int j = 0; j < nMode_; ++j) {
               totalFourierModes_(i, j) = fourierModes_(i, j);
            }
         }
         #endif

         if (simulation().domain().isMaster()) {
            // Increment structure factors
            double volume = simulation().boundary().volume();
            double norm;
            for (j = 0; j < nMode_; ++j) {
               double maxValue = 0.0;
               IntVector maxIntVector;
               double maxQ;
               for (i = 0; i < nWave_; ++i) {
                  norm = std::norm(totalFourierModes_(i, j));
                  if ( double(norm/volume) >= maxValue ) {
                     maxValue = double(norm/volume);
                     maxIntVector = waveIntVectors_[i];
                     maxQ = waveVectors_[i].abs();
                  }
                  structureFactors_(i, j) += norm/volume;
               }
               maximumValue_[j].insert(maximumValue_[j].end(), 1, maxValue);
               maximumWaveIntVector_[j].insert(maximumWaveIntVector_[j].end(), 
                                               1, maxIntVector);
               maximumQ_[j].insert(maximumQ_[j].end(), 1, maxQ);
            }
         }

         ++nSample_;
      }

   }

   /*
   * Calculate floating point wavevectors.
   */
   void StructureFactor::makeWaveVectors() 
   {
      Vector    dWave;
      Boundary* boundaryPtr = &simulation().boundary();
      int       i, j;

      // Calculate wavevectors
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
   * Write data to three output files.
   */
   void StructureFactor::output()
   {
      if (simulation().domain().isMaster()) {

         // Write parameters to a *.prm file
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), 
                                                  outputFile_);
         writeParam(outputFile_);
         outputFile_.close();

         // Output structure factors to one *.dat file
         simulation().fileMaster().openOutputFile(outputFileName(".dat"), 
                                                  outputFile_);
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

         // Output maximum structure factors to *_max.dat file
         simulation().fileMaster().openOutputFile(outputFileName("_max.dat"), 
                                                  outputFile_);
         for (j = 0; j < nMode_; ++j) {
            for (int i = 0; i < nSample_; ++i) {
               outputFile_ << maximumWaveIntVector_[j][i];
               outputFile_ << Dbl(maximumQ_[j][i], 20, 8);
               outputFile_ << Dbl(maximumValue_[j][i], 20, 8);
               outputFile_ << std::endl;
            }
         }
         outputFile_.close();
      }

   } 

}
#endif
