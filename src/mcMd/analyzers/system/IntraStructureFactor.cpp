#ifndef MCMD_INTRA_STRUCTURE_FACTOR_CPP
#define MCMD_INTRA_STRUCTURE_FACTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraStructureFactor.h"
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
   IntraStructureFactor::IntraStructureFactor(System& system) 
    : SystemAnalyzer<System>(system),
      isInitialized_(false)
   {  setClassName("IntraStructureFactor"); }

   /*
   * Destructor.
   */
   IntraStructureFactor::~IntraStructureFactor()
   {}

   /*
   * Read parameters from file, and allocate data array.
   */
   void IntraStructureFactor::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);

      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "nAtomTypeIdPair", nAtomTypeIdPair_);
      atomTypeIdPairs_.allocate(nAtomTypeIdPair_);
      readDArray< Pair<int> >(in, "atomTypeIdPairs", atomTypeIdPairs_, 
                                                     nAtomTypeIdPair_);
      read<int>(in, "nWave", nWave_);

      waveIntVectors_.allocate(nWave_);
      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, nWave_);
      nAtomType_ = system().simulation().nAtomType();

      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nAtomType_ + 1);
      structureFactors_.allocate(nWave_, nAtomTypeIdPair_);
      structureFactorDelta_.allocate(nWave_, nAtomTypeIdPair_);

      isInitialized_ = true; 
   }

   /*
   * Load state from an archive.
   */
   void IntraStructureFactor::loadParameters(Serializable::IArchive& ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "nAtomTypeIdPair", nAtomTypeIdPair_);
      atomTypeIdPairs_.allocate(nAtomTypeIdPair_);
      loadDArray< Pair<int> >(ar, "atomTypeIdPairs", atomTypeIdPairs_, 
                                                     nAtomTypeIdPair_);
      loadParameter<int>(ar, "nWave", nWave_);

      //waveIntVectors_.allocate(nWave_);
      loadDArray<IntVector>(ar, "waveIntVectors", waveIntVectors_, nWave_);

      ar & nAtomType_;
      if (nAtomType_ != system().simulation().nAtomType()) {
         UTIL_THROW("Inconsistent nAtomType");
      }
      // structureFactors_.allocate(nWave_, nAtomTypeIdPair_);
      ar & structureFactors_;
      //fourierModes_.allocate(nWave_, nAtomType_ + 1);
      ar & fourierModes_;
      ar & nSample_;

      waveVectors_.allocate(nWave_);
      isInitialized_ = true; 
   }


   /*
   * Save state to archive.
   */
   void IntraStructureFactor::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Set up before simulation.
   */
   void IntraStructureFactor::setup() 
   {
      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("Object is not initialized");
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

   /* 
   * Increment Structure Factor
   */
   void IntraStructureFactor::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         Vector  position;
         std::complex<double> rho[2];
         std::complex<double> expFactor;
         double  volume = system().boundary().volume();
         double  product, dS;
         System::ConstMoleculeIterator molIter;
         Molecule::ConstAtomIterator   atomIter;
         int  typeId, i, j, k;

         makeWaveVectors();

         // Set all Deltas to zero
         for (i = 0; i < nWave_; ++i) {
            for (int pairId = 0; pairId < nAtomTypeIdPair_; ++pairId) {
               structureFactorDelta_(i, pairId) = 0;
            }
         }
 
         // Loop over molecules in species
         system().begin(speciesId_, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {

            // Set all Fourier modes to zero
            for (i = 0; i < nWave_; ++i) {
               for (typeId = 0; typeId <= nAtomType_; ++typeId) {
                  fourierModes_(i, typeId) = 0;
               }
            }
 
            // Loop over all atoms in one molecule
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

            // Increment structure factors
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
                  dS = std::real(rho[0]*rho[1])/volume;
                  structureFactors_(i, j) += dS;
                  structureFactorDelta_(i,j) += dS;
               }
            }

         }

         ++nSample_;
      }

   }

   /**
   * Calculate floating point wavevectors.
   */
   void IntraStructureFactor::makeWaveVectors() 
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

   void IntraStructureFactor::output() 
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

   }

}
#endif
