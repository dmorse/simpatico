#ifndef MCMD_TYPE_DISTRIBUTION_CPP
#define MCMD_TYPE_DISTRIBUTION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TypeDistribution.h"

#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/species/SpeciesMutator.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <mcMd/util/FileMaster.h>

namespace McMd
{

   using namespace Util;

   // Constructor
   TypeDistribution::TypeDistribution(McSystem& system)
    : SystemDiagnostic<McSystem>(system),
      outputFile_(),
      distribution_(),
      nSample_(0),
      speciesId_(-1),
      nState_(-1),
      speciesPtr_(0),
      mutatorPtr_(0)
   {  setClassName("TypeDistribution"); }

   void TypeDistribution::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);

      speciesPtr_ = &(system().simulation().species(speciesId_));
      if (!speciesPtr_->isMutable()) {
         UTIL_THROW("Error: Species must be mutable");
      }
      mutatorPtr_ = &speciesPtr_->mutator();
      nState_ = mutatorPtr_->nState();

      // Allocate IntDistribution accumulator
      distribution_.allocate(nState_);
      for (int i = 0; i < nState_; ++i) {
         distribution_[i] = 0;
      }

   }
 
   // Evaluate energy and print.
   void TypeDistribution::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {

         for (int i = 0; i < nState_; ++i) {
            distribution_[i] += mutatorPtr_->stateOccupancy(i);
         }
         ++nSample_;

      }
   }
 
   /*
   * Summary
   */
   void TypeDistribution::output() 
   {
      // Open and write summary file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);


      double norm = 1.0/double(nSample_);
      for (int i = 0; i < nState_; ++i) {
         outputFile_ << Int(i,5) 
                     << Dbl(distribution_[i]*norm, 15) << std::endl;
      }
      outputFile_.close();
   }
   
}
#endif 
