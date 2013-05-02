#ifndef MCMD_VAN_HOVE_CPP
#define MCMD_VAN_HOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "VanHove.h"
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
   VanHove::VanHove(System& system) 
    : SystemDiagnostic<System>(system),
      isInitialized_(false)
   {  setClassName("VanHove"); }

   /*
   * Destructor
   */
   VanHove::~VanHove() 
   {}

   /*
   * Read parameters from file, and allocate memory.
   */
   void VanHove::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);

      nAtomType_ = system().simulation().nAtomType();
      atomTypeCoeffs_.allocate(nAtomType_);
      readDArray<double>(in, "atomTypeCoeffs", atomTypeCoeffs_, nAtomType_);

      read<int>(in, "nBuffer", nBuffer_);
      read<int>(in, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_);
      accumulators_.allocate(nWave_);
      for (int i = 0; i < nWave_; ++i) {
         accumulators_[i].setParam(nBuffer_);
      }

      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, nWave_);
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void VanHove::loadParameters(Serializable::IArchive& ar)
   {  
      Diagnostic::loadParameters(ar);
      ar & nAtomType_;
      atomTypeCoeffs_.allocate(nAtomType_);
      loadDArray<double>(ar, "atomTypeCoeffs", atomTypeCoeffs_, nAtomType_);
      loadParameter<int>(ar, "nBuffer", nBuffer_);
      loadParameter<int>(ar, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      loadDArray<IntVector>(ar, "waveIntVectors", waveIntVectors_, nWave_);
      accumulators_.allocate(nWave_);
      ar & accumulators_;
      ar & nSample_;

      // Validate
      if (nAtomType_ != system().simulation().nAtomType()) {
         UTIL_THROW("Inconsistent values of nAtomType");
      }
      if (nAtomType_ != atomTypeCoeffs_.capacity()) {
         UTIL_THROW("Inconsistent capacity for atomTypeCoeffs");
      }
      if (nWave_ != waveIntVectors_.capacity()) {
         UTIL_THROW("Inconsistent capacity for waveIntVectors");
      }

      // Allocate work arrays.
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_);

      isInitialized_ = true;
   }

   /*
   * Save state to an archive.
   */
   void VanHove::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /*
   * Clear accumulators.
   */
   void VanHove::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized");
      }
      assert (nBuffer_ > 0);
      assert (nWave_ > 0);
      for (int i = 0; i < nWave_; ++i) {
         accumulators_[i].clear();
      }
   }
 
   /// Increment Structure Factor
   void VanHove::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         Vector  position;
         std::complex<double>  expFactor;
         double  product, coeff;
         System::ConstMoleculeIterator  molIter;
         Molecule::ConstAtomIterator  atomIter;
         int  nSpecies, iSpecies, i;

         makeWaveVectors();

         // Set all Fourier modes to zero
         for (i = 0; i < nWave_; ++i) {
            fourierModes_[i] = std::complex<double>(0.0, 0.0);
         }
 
         // Loop over all atoms to calculate Fourier modes
         nSpecies = system().simulation().nSpecies();
         for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            system().begin(iSpecies, molIter); 
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(atomIter); 
               for ( ; atomIter.notEnd(); ++atomIter) {
                  position = atomIter->position();
                  coeff    = atomTypeCoeffs_[atomIter->typeId()];
 
		  // Loop over wavevectors
                  for (i = 0; i < nWave_; ++i) {
               
                     product = position.dot(waveVectors_[i]);
                     expFactor = exp( product*Constants::Im );
                     fourierModes_[i] += (coeff*expFactor);
		 
                  }
  
               }
            }

         }

         // Add Fourier modes to autocorrelation accumulators
         double sqrtV = sqrt(system().boundary().volume());
         for (i = 0; i < nWave_; ++i) {
            accumulators_[i].sample(fourierModes_[i]/sqrtV);
         }

      }

   }

   /**
   * Calculate floating point wavevectors.
   */
   void VanHove::makeWaveVectors() 
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
         }
      }
   }

   void VanHove::output() 
   {
      std::string suffix;
      int i, j;

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      // Output autocorrelation functions to file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      for (i = 0; i < nWave_; ++i) {

         #if 0
         // Construct file suffix for this structure factor
         suffix = std::string("_");
         for (j = 0; j < Dimension; ++j) {
            suffix += toString(waveIntVectors_[i][j]);
         }
         suffix += std::string(".dat");
         outputFile_ << suffix << std::endl;
         //fileMaster().openOutputFile(outputFileName(suffix), outputFile_);
         #endif

         for (j = 0; j < Dimension; ++j) {
            outputFile_ << Int(waveIntVectors_[i][j], 5);
         }
         outputFile_ << Dbl(waveVectors_[i].abs(), 20, 8);
         outputFile_ << std::endl;
         accumulators_[i].output(outputFile_);
         outputFile_ << std::endl;
         outputFile_ << std::endl;
         //outputFile_.close();

      }
      outputFile_.close();

   }

}
#endif
