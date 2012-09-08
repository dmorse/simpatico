#ifndef MCMD_COMPOSITION_PROFILE_CPP
#define MCMD_COMPOSITION_PROFILE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CompositionProfile.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/McMd_mpi.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Dimension.h>
#include <mcMd/util/FileMaster.h>
#include <util/util/ioUtil.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   CompositionProfile::CompositionProfile(System& system) 
    : SystemDiagnostic<System>(system),
      isInitialized_(false)
   {  setClassName("CompositionProfile"); }

   CompositionProfile::~CompositionProfile() 
   {}

   /// Read parameters from file, and allocate direction vectors.
   void CompositionProfile::readParameters(std::istream& in) 
   {

      // Read interval and output file 
      readInterval(in);
      readOutputFileName(in);

      nAtomType_ = system().simulation().nAtomType();

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
   void CompositionProfile::setup() 
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
      Vector blengths;

      // Lengths of simulation box
      blengths = system().boundary().lengths();

      if (isAtInterval(iStep))  {

         Vector  position;
         double  product;
         System::ConstMoleculeIterator  molIter;
         Molecule::ConstAtomIterator  atomIter;
         int  nSpecies, iSpecies, typeId;

         // Loop over all atoms
         nSpecies = system().simulation().nSpecies();
         for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            system().begin(iSpecies, molIter); 
          
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(atomIter); 
          
               for ( ; atomIter.notEnd(); ++atomIter) {
                  position = atomIter->position();
                  typeId   = atomIter->typeId();
                  for (int i = 0; i < Dimension; ++i) {
                     position[i] /= blengths[i];
                  }
		  // Loop over direction vectors
                  for (int i = 0; i < nDirection_; ++i) {
                     product = position.dot(waveVectors_[i]);
                     product /= waveVectors_[i].abs();
                     accumulators_[i+typeId*nDirection_].sample(product);
                  }
		 
               }
  
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
      Boundary* boundaryPtr = &system().boundary();
      int       i, j;

      // Calculate wavevectors
      for (i = 0; i < nDirection_; ++i) {
         waveVectors_[i] = Vector::Zero;
         for (j = 0; j < Dimension; ++j) {
            dWave  = boundaryPtr->reciprocalBasisVector(j);
            dWave *= intVectors_[i][j];
            waveVectors_[i] += dWave;
            //std::cout << waveVectors_[i] << std::endl;
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
