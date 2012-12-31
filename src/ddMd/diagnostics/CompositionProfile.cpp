#ifndef DDMD_COMPOSITION_PROFILE_CPP
#define DDMD_COMPOSITION_PROFILE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CompositionProfile.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/simulation/SimulationAccess.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <util/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Dbl.h>

namespace DdMd
{

   using namespace Util;

   /// Constructor.
   CompositionProfile::CompositionProfile(Simulation& simulation) 
    : Diagnostic,
      isInitialized_(false)
   { }

   CompositionProfile::~CompositionProfile() 
   {}

   /// Read parameters from file, and allocate direction vectors.
   void CompositionProfile::readParameters(std::istream& in) 
   {

      // Read interval and output file 
      readInterval(in);
      readOutputFileName(in);

      nAtomType_ = simulation().nAtomType();

      // Read number of direction vectors and direction vectors 
      read<int>(in, "nDirections", nDirection_);
      intVectors_.allocate(nDirection_);
      waveVectors_.allocate(nDirection_);
      readDArray<IntVector>(in, "intVectors", intVectors_, nDirection_);

      accumulators_.allocate(nDirection_*nAtomType_);

      isInitialized_ = true;
   
   }

   /*
   * Clear accumulators.
   */
   void CompositionProfile::clear() 
   {  
      int bin;
      double min, max;
      
      for (int i=0; i < nDirection_; ++i) {
         for (int j=0; j < nAtomType_; ++j) {

            min = 0.0;
            max = 1.0;

            // number of bins in histogram = int(((max-min)/0.05));
            bin = 500;

            accumulators_[i+j*nDirection_].setParam(min, max, bin);
         }
      }

      makeWaveVectors();

      // Clear accumulators
      for (int i = 0; i < nDirection_; ++i) {
         for (int j = 0; j < nAtomType_; ++j){
            accumulators_[i+j*nDirection_].clear();
         }
      } 
      nSample_ = 0;
   }
 
   void CompositionProfile::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         
         Vector blengths;

         // Lengths of simulation box
         blengths = simulation().boundary().lengths();

         Vector  position;
         double  product;
         System::ConstMoleculeIterator  molIter;
         Molecule::ConstAtomIterator  atomIter;
         int  nSpecies, iSpecies, typeId;

         // Loop over all atoms
         simulation().atomStorage().begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            position = atomIter->position();
            typeId   = atomIter->typeId();

            for (int i = 0; i < Dimension; ++i) {
               position[i] /= blengths[i];
            }

            for (int i = 0; i < nDirection_; ++i) {
               product = position.dot(waveVectors_[i]);
               product /= waveVectors_[i].abs();
               accumulators_[i+typeId*nDirection_].sample(product);
            }
         }

         if (simulation().domain().isMaster()) {
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
               maximumWaveIntVector_[j].insert(maximumWaveIntVector_[j].end(), 1, maxIntVector);
               maximumQ_[j].insert(maximumQ_[j].end(), 1, maxQ);
            }
         }


         ++nSample_;
      }
   }

   /**
   * Calculate floating point wavevectors.
   */
   void CompositionProfile::makeWaveVectors()
   {
      Vector    dWave;
      Boundary* boundaryPtr = &simulation().boundary();
      int       i, j;

      // Calculate wavevectors
      for (i = 0; i < nDirection_; ++i) {
         waveVectors_[i] = Vector::Zero;
         for (j = 0; j < Dimension; ++j) {
            dWave  = boundaryPtr->reciprocalBasisVector(j);
            dWave *= intVectors_[i][j];
            waveVectors_[i] += dWave;
         }
      }
   }

   void CompositionProfile::output() 
   {
      int i, j, k;
      std::string suffix;
      
      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      // Output statistical analysis to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      for (i = 0; i < nDirection_; ++i) {
         for (j = 0; j < Dimension; ++j) {
            outputFile_ << Dbl(waveVectors_[i][j], 5);
         }
         outputFile_ << std::endl;
      
         for (k = 0; k < nAtomType_; ++k) {
            accumulators_[i+k*nDirection_].output(outputFile_);
            outputFile_ << std::endl;
            outputFile_ << std::endl;
         }
      }
      outputFile_.close();


   }

}
#endif
