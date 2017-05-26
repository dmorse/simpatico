/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SemiGrandDistribution.h"

#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/SpeciesMutator.h>
#include <simp/species/Species.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/misc/FileMaster.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   SemiGrandDistribution::SemiGrandDistribution(McSystem& system)
    : SystemAnalyzer<McSystem>(system),
      outputFile_(),
      distribution_(),
      speciesId_(-1),
      moleculeCapacity_(-1),
      speciesPtr_(0),
      mutatorPtr_(0)
   {  setClassName("SemiGrandDistribution"); }

   /*
   * Read parameters from file and initialize. 
   */
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

   /*
   * Load state from an archive.
   */
   void SemiGrandDistribution::loadParameters(Serializable::IArchive& ar)
   {  
      Analyzer::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      ar & moleculeCapacity_;
      ar & distribution_;

      // Validate values and set pointer members.
      if (speciesId_ < 0 || speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("Invalid speciesId");
      }
      if (moleculeCapacity_ <= 0) {
         UTIL_THROW("Error: moleculeCapacity <= 0");
      }
      speciesPtr_ = &(system().simulation().species(speciesId_));
      if (moleculeCapacity_ != speciesPtr_->capacity()) {
          UTIL_THROW("Inconsiste moleculeCapacity");
      }
      if (!speciesPtr_->isMutable()) {
         UTIL_THROW("Error: Species must be mutable");
      }
      mutatorPtr_ = &speciesPtr_->mutator();
   }

   /*
   * Save state to an archive.
   */
   void SemiGrandDistribution::save(Serializable::OArchive& ar)
   {  ar & *this; }


   /* 
   * Sample occupation.
   */
   void SemiGrandDistribution::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {
         distribution_.sample(mutatorPtr_->stateOccupancy(0));
      }
   }
 
   /*
   * Open and write summary output file.
   */
   void SemiGrandDistribution::output() 
   {
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      distribution_.output(outputFile_);
      outputFile_.close();
   }
   
}
