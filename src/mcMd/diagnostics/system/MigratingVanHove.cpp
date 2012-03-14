#ifdef UTIL_MPI
#ifdef MCMD_PERTURB

#ifndef MIGRATING_VAN_HOVE_CPP
#define MIGRATING_VAN_HOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MigratingVanHove.h"
#include <mcMd/perturb/Perturbation.h>
#include <mcMd/perturb/LinearPerturbation.h>
#include <mcMd/perturb/ReplicaMove.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/McMd_mpi.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/perturb/ReplicaMove.h>
#include <util/math/Constants.h>
#include <util/space/Dimension.h>
#include <mcMd/util/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/util/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   MigratingVanHove::MigratingVanHove(System& system) 
    : SystemDiagnostic<System>(system),
      communicatorPtr_(0),
      replicaMovePtr_(&system.replicaMove()),
      isInitialized_(false)
   {  
      communicatorPtr_ = &(system.simulation().communicator());
   }
  
 
   MigratingVanHove::~MigratingVanHove() 
   {}

   /// Read parameters from file, and allocate data array.
   void MigratingVanHove::readParam(std::istream& in) 
   {

      // Read interval and parameters for AutoCorrArray
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
   * Clear accumulators.
   */
   void MigratingVanHove::initialize() 
   {
      replicaMovePtr_->registerObserver(*this);
      assert (nBuffer_ > 0);
      assert (nWave_ > 0);
      for (int i = 0; i < nWave_; ++i) {
         accumulators_[i].clear();
      }
   }
 
   /// Increment Structure Factor
   void MigratingVanHove::sample(long iStep) 
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
   void MigratingVanHove::makeWaveVectors() 
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

   void MigratingVanHove::output() 
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

   void MigratingVanHove::update(const int &partnerId)
   {
      int size;
      size = memorySize(accumulators_);
      size += memorySize(fourierModes_);
      size += memorySize(nBuffer_);
      size += memorySize(nSample_);
      
      MemoryOArchive sendCurrent;
      sendCurrent.allocate(size);
      
      sendCurrent << accumulators_;
      sendCurrent << fourierModes_;
      sendCurrent << nBuffer_;
      sendCurrent << nSample_;
      
      sendCurrent.send(*communicatorPtr_, partnerId);
       
      MemoryIArchive recvCurrent;
      recvCurrent.allocate(size);
 
      recvCurrent.recv(*communicatorPtr_, partnerId);
      
      recvCurrent >> accumulators_;
      recvCurrent >> fourierModes_;
      recvCurrent >> nBuffer_;
      recvCurrent >> nSample_;
   }

   /*
   * Save state to binary file archive.
   */
   void MigratingVanHove::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void MigratingVanHove::load(Serializable::IArchiveType& ar)
   { ar & *this; }


}
#endif

#endif // UTIL_MPI
#endif
