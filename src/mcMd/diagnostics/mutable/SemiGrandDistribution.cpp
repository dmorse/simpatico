#ifndef MCMD_SEMI_GRAND_DISTRIBUTION_CPP
#define MCMD_SEMI_GRAND_DISTRIBUTION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SemiGrandDistribution.h"

#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/species/SpeciesMutator.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <mcMd/misc/FileMaster.h>

namespace McMd
{

   using namespace Util;

   // Constructor
   SemiGrandDistribution::SemiGrandDistribution(McSystem& system)
    : SystemDiagnostic<McSystem>(system),
      outputFile_(),
      distribution_(),
      speciesId_(-1),
      moleculeCapacity_(-1),
      speciesPtr_(0),
      mutatorPtr_(0)
   { setClassName("SemiGrandDistribution"); }

   void SemiGrandDistribution::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);

      speciesPtr_ = &(system().simulation().species(speciesId_));
      moleculeCapacity_ = speciesPtr_->capacity();
      if (moleculeCapacity_ <= 0) {
         UTIL_THROW("Error: moleculeCapacity <= 0");
      }
      if (!speciesPtr_->isMutable()) {
         UTIL_THROW("Error: Species must be mutable");
      }
      mutatorPtr_ = &speciesPtr_->mutator();
     
      // Allocate IntDistribution accumulator
      distribution_.setParam(0, moleculeCapacity_);
      distribution_.clear();

   }
 
   // Sample occupation
   void SemiGrandDistribution::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {
         distribution_.sample(mutatorPtr_->stateOccupancy(0));
      }
   }
 
   /*
   * Summary
   */
   void SemiGrandDistribution::output() 
   {
      // Open and write summary file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      distribution_.output(outputFile_);
      //double norm = 1.0/double(nSample_);
      //for (int i = 0; i < moleculeCapacity_; ++i) {
      //   outputFile_ << Int(i,5) 
      //               << Dbl(distribution_[i]*norm, 15) << std::endl;
      //}
      outputFile_.close();
   }
   
}
#endif 
