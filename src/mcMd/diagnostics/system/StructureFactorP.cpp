#ifndef MCMD_STRUCTURE_FACTOR_P_CPP
#define MCMD_STRUCTURE_FACTOR_P_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "StructureFactorP.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/McMd_mpi.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/math/Constants.h>
#include <util/space/Dimension.h>
#include <mcMd/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   StructureFactorP::StructureFactorP(System& system) 
    : SystemDiagnostic<System>(system)
   {  setClassName("StructureFactorP"); }

   StructureFactorP::~StructureFactorP() 
   {}

   /// Read parameters from file, and allocate data array.
   void StructureFactorP::readParameters(std::istream& in) 
   {

      // Read interval and parameters for AutoCorrArray
      //SystemDiagnostic<System>::readParameters(in);
      readInterval(in);
      readOutputFileName(in);

      nAtomType_ = system().simulation().nAtomType();

      read<int>(in, "nAtomTypeIdPair", nAtomTypeIdPair_);
      atomTypeIdPairs_.allocate(nAtomTypeIdPair_);
      readDArray< Pair<int> >(in, "atomTypeIdPairs", atomTypeIdPairs_, 
                                                     nAtomTypeIdPair_);

      read<int>(in, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, nWave_);

      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nAtomType_ + 1);
      structureFactors_.allocate(nWave_, nAtomTypeIdPair_);

      isInitialized_ = true;
   }

   /*
   * Setup immediately before simulation.
   */
   void StructureFactorP::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Object not initialized");
      }
      assert (nWave_ > 0);
      assert (nAtomTypeIdPair_ > 0);

      // Clear accumulators
      int i, j;
      for (i = 0; i < nWave_; ++i) {
         for (j = 0; j < nAtomTypeIdPair_; ++j) {
            structureFactors_(i, j) = 0.0;
         }
      }
      nSample_ = 0;

   }
 
   /// Increment Structure Factor
   void StructureFactorP::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         Vector  position;
         std::complex<double> expFactor;
         double  product;
         System::ConstMoleculeIterator  molIter;
         Molecule::ConstAtomIterator  atomIter;
         int  nSpecies, iSpecies, typeId, i;

         makeWaveVectors();

         // Set all Fourier modes to zero
         for (i = 0; i < nWave_; ++i) {
            for (typeId = 0; typeId <= nAtomType_; ++typeId) {
               fourierModes_(i, typeId) = std::complex<double>(0.0, 0.0);
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
                     fourierModes_(i, typeId)     += expFactor;
                     fourierModes_(i, nAtomType_) += expFactor;
		 
                  }
  
               }
            }

         }

         // Increment structure factors
         std::complex<double> rho[2];
         double  volume = system().boundary().volume();
         int     j, k;
         for (i = 0; i < nWave_; ++i) {
            for (j = 0; j < nAtomTypeIdPair_; ++j) {
               for (k = 0; k < 2; ++k) {
                  typeId = atomTypeIdPairs_[j][k];
                  if (typeId >= 0) {
                     rho[k] = fourierModes_(i, typeId);
                  } else {
                     rho[k] = fourierModes_(i, nAtomType_);
                  }
               }
               rho[0] = std::conj(rho[0]);
               structureFactors_(i, j) += std::real(rho[0]*rho[1])/volume;
            }
         }

         ++nSample_;

      }

   }

   /**
   * Calculate floating point wavevectors.
   */
   void StructureFactorP::makeWaveVectors() 
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

   void StructureFactorP::output() 
   {
      double      value;
      int         i, j, k;

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output structure factors to one file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      for (i = 0; i < nWave_; ++i) {
         for (k = 0; k < Dimension; ++k) {
            outputFile_ << Int(waveIntVectors_[i][k], 5);
         }
         outputFile_ << Dbl(waveVectors_[i].abs(), 20, 8);
         for (j = 0; j < nAtomTypeIdPair_; ++j) {
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
      for (j = 0; j < nAtomTypeIdPair_; ++j) {

         // Construct file suffix for this structure factor
         suffix = std::string(".");
         for (k = 0; k < 2; k++) {
            typeId = atomTypeIdPairs_[j][k];
            if (typeId < 0) {
               suffix += std::string("A");
            } else {
               suffix += toString(typeId);
            }
         }
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
