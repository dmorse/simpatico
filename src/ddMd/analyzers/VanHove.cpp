#ifndef DDMD_VAN_HOVE_CPP
#define DDMD_VAN_HOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "VanHove.h"
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
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   VanHove::VanHove(Simulation& simulation) 
    : Analyzer(simulation),
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

      nAtomType_ = simulation().nAtomType();
      atomTypeCoeffs_.allocate(nAtomType_);
      readDArray<double>(in, "atomTypeCoeffs", atomTypeCoeffs_, nAtomType_);

      read<int>(in, "nBuffer", nBuffer_);
      read<int>(in, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_);
      totalFourierModes_.allocate(nWave_);

      if (simulation().domain().isMaster()) {
         accumulators_.allocate(nWave_);
         for (int i = 0; i < nWave_; ++i) {
            accumulators_[i].setParam(nBuffer_);
         }
      }

      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, nWave_);
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void VanHove::loadParameters(Serializable::IArchive& ar)
   {  
      nAtomType_ = simulation().nAtomType();
      loadInterval(ar);
      loadOutputFileName(ar);
      atomTypeCoeffs_.allocate(nAtomType_);
      loadDArray<double>(ar, "atomTypeCoeffs", atomTypeCoeffs_, nAtomType_);
      loadParameter<int>(ar, "nBuffer", nBuffer_);
      loadParameter<int>(ar, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      loadDArray<IntVector>(ar, "waveIntVectors", waveIntVectors_, nWave_);

      // Load and broadcast nSample_
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);
      
      if (simulation().domain().isMaster()) {
         accumulators_.allocate(nWave_);
         ar >> accumulators_;
      }
             
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_);
      totalFourierModes_.allocate(nWave_);

      isInitialized_ = true;
   }

   /*
   * Save state to an archive.
   */
   void VanHove::save(Serializable::OArchive& ar)
   { 
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << atomTypeCoeffs_;
      ar << nBuffer_;
      ar << nWave_;
      ar << waveIntVectors_;

      ar << nSample_;

      if (simulation().domain().isMaster()) {
         ar << accumulators_;
      }
   }

   /*
   * Clear accumulators.
   */
   void VanHove::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized");
      }
      assert (nBuffer_ > 0);
      assert (nWave_ > 0);
      nSample_ = 0;
      if (simulation().domain().isMaster()) {
         for (int i = 0; i < nWave_; ++i) {
            accumulators_[i].clear();
         }
      }
   }
 
   /// Increment Structure Factor
   void VanHove::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         Vector  position;
         std::complex<double>  expFactor;
         double  product, coeff;
         AtomIterator  atomIter;
         int  i, typeId;

         makeWaveVectors();

         // Set all Fourier modes to zero
         for (i = 0; i < nWave_; ++i) {
            fourierModes_[i] = std::complex<double>(0.0, 0.0);
         }
 
         simulation().atomStorage().begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            position = atomIter->position();
            typeId   = atomIter->typeId();
            coeff    = atomTypeCoeffs_[typeId];
            
            // Loop over wavevectors
            for (i = 0; i < nWave_; ++i) {

               product = position.dot(waveVectors_[i]);
               expFactor = exp( product*Constants::Im );
               fourierModes_[i] += coeff*expFactor;
            }
         }
  
         for (i = 0; i < nWave_; ++i) {
            totalFourierModes_[i] = std::complex<double>(0.0, 0.0);
         }

         #ifdef UTIL_MPI
         // Loop over wavevectors
         for (int i = 0; i < nWave_; ++i) {
            //Sum values from all processors.
            simulation().domain().communicator().
            Reduce(&fourierModes_[i], &totalFourierModes_[i],
                    1, MPI::DOUBLE_COMPLEX, MPI::SUM, 0);
         }
         #else
         for (int i = 0; i < nWave_; ++i) {
            totalFourierModes_[i] = fourierModes_[i];
         }
         #endif
         
         if (simulation().domain().isMaster()) {
            // Add Fourier modes to autocorrelation accumulators
            double sqrtV = sqrt(simulation().boundary().volume());
            for (int i = 0; i < nWave_; ++i) {
               accumulators_[i].sample(totalFourierModes_[i]/sqrtV);
            }
         }
         ++nSample_;
      }

   }

   /**
   * Calculate floating point wavevectors.
   */
   void VanHove::makeWaveVectors() 
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

   void VanHove::output() 
   {
      if (simulation().domain().isMaster()) {
         // Write parameters to a *.prm file
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         writeParam(outputFile_);
         outputFile_.close();
        
         // Output structure factors to one *.dat file
         simulation().fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
        
         int i, j;

         // Output autocorrelation functions to file
         for (i = 0; i < nWave_; ++i) {

            for (j = 0; j < Dimension; ++j) {
               outputFile_ << Int(waveIntVectors_[i][j], 5);
            }
            outputFile_ << Dbl(waveVectors_[i].abs(), 20, 8);
            outputFile_ << std::endl;
            accumulators_[i].output(outputFile_);
            outputFile_ << std::endl;
            outputFile_ << std::endl;
         }
         outputFile_.close();
      }

   }

}
#endif
