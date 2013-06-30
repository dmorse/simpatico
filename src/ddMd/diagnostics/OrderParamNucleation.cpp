#ifndef DDMD_ORDER_PARAM_NUCLEATION_CPP
#define DDMD_ORDER_PARAM_NUCLEATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OrderParamNucleation.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/simulation/SimulationAccess.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <util/boundary/Boundary.h>
#include <util/math/Constants.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace DdMd
{

   using namespace Util;

   /// Constructor.
   OrderParamNucleation::OrderParamNucleation(Simulation& simulation) 
    : Diagnostic(simulation),
      isInitialized_(false)
   {  setClassName("OrderParamNucleation"); }

   OrderParamNucleation::~OrderParamNucleation() 
   {}

   /// Read parameters from file, and allocate data array.
   void OrderParamNucleation::readParameters(std::istream& in) 
   {
      nAtomType_ = simulation().nAtomType();

      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "perpIndex", perpIndex_);
      if (perpIndex_ < 0 || perpIndex_ >= Dimension) {
         UTIL_THROW("Invalid index for ordering direction.");
      }
      read<int>(in, "parallelIndex", parallelIndex_);
      if (parallelIndex_ < 0 || parallelIndex_ >= Dimension) {
         UTIL_THROW("Invalid index for parallel direction.");
      }
      read<int>(in, "periodicity", periodicity_);
      read<int>(in, "nBin", nBin_);

      cosFactors_.allocate(nAtomType_, nBin_);
      totalCosFactors_.allocate(nAtomType_, nBin_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void OrderParamNucleation::loadParameters(Serializable::IArchive &ar)
   {
      nAtomType_ = simulation().nAtomType();

      // Load and broadcast parameter file parameters
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar, "perpIndex", perpIndex_);
      if (perpIndex_ < 0 || perpIndex_ >= Dimension) {
         UTIL_THROW("Invalid index for ordering direction.");
      }
      loadParameter<int>(ar, "parallelIndex", parallelIndex_);
      if (parallelIndex_ < 0 || parallelIndex_ >= Dimension) {
         UTIL_THROW("Invalid index for parallel direction.");
      }
      loadParameter<int>(ar, "periodicity", periodicity_);
      loadParameter<int>(ar, "nBin", nBin_);

      cosFactors_.allocate(nAtomType_, nBin_);
      totalCosFactors_.allocate(nAtomType_, nBin_);

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void OrderParamNucleation::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << perpIndex_;
      ar << parallelIndex_;
      ar << periodicity_;
      ar << nBin_;
   }
  
   /*
   * Clear accumulators.
   */
   void OrderParamNucleation::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
      nSample_ = 0;

      int i, j;
      for (i = 0; i < nAtomType_; ++i) {
         for (j = 0; j < nBin_; ++j) {
            cosFactors_(i, j) = 0.0;
            totalCosFactors_(i, j) = 0.0;
         }
      }

   }

   /*
   * Increment cos factor.
   */
   void OrderParamNucleation::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         Vector blengths;
         Boundary* boundaryPtr = &simulation().boundary();
         blengths = boundaryPtr->lengths();

         Vector position;
         double  slabLength, product, cosFactor;
         AtomIterator  atomIter;
         int typeId, bin;

         simulation().atomStorage().begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            position = atomIter->position();
            typeId   = atomIter->typeId();
 
            slabLength = blengths[parallelIndex_]/nBin_;
            bin = position[parallelIndex_]/slabLength;

            product = position[perpIndex_]*2.0*Constants::Pi*periodicity_;
            product /= blengths[perpIndex_];
            cosFactor = cos(product)*cos(product);
            cosFactors_(typeId, bin) += cosFactor;
         }

         #ifdef UTIL_MPI
         // Loop over wavevectors
         for (int i = 0; i < nAtomType_; ++i) {
            for (int j = 0; j < nBin_; ++j) {
            //Sum values from all processors.
            simulation().domain().communicator().
                         Reduce(&cosFactors_(i, j), &totalCosFactors_(i, j),
                                1, MPI::DOUBLE, MPI::SUM, 0);
            }
         }
         #else
         for (int i = 0; i < nAtomType_; ++i) {
            for (int j = 0; j < nBin_; ++j) {
               totalCosFactors_(i, j)  = cosFactors_(i, j);
            }
         }
         #endif

         ++nSample_;
      }

   }

   /*
   * Write data to three output files.
   */
   void OrderParamNucleation::output()
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
         double      volume, value;
         int         i, j;
         volume = simulation().boundary().volume();
         for (i = 0; i < nBin_; ++i) {
            outputFile_ << Int(i, 4);
            for (j = 0; j < nAtomType_; ++j) {
               value = totalCosFactors_(j, i)/volume/double(nSample_);
               outputFile_ << Dbl(value, 18, 8);
            }
            outputFile_ << std::endl;
         }
         outputFile_.close();
      }

   } 

}
#endif
