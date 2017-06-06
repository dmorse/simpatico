/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BondLengthDist.h"
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/math/feq.h>
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   BondLengthDist::BondLengthDist(System& system) 
    : SystemAnalyzer<System>(system),
      min_(0.0),
      max_(0.0),
      nBin_(0),
      isInitialized_(false)
   {  setClassName("BondLengthDist"); }

   /// Read parameters from file, and allocate data array.
   void BondLengthDist::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      read<double>(in, "min", min_);
      read<double>(in, "max", max_);
      read<double>(in, "nBin", nBin_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      //readParamComposite(in, accumulator_);
      accumulator_.setParam(min_, max_, nBin_);
      accumulator_.clear();
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void BondLengthDist::loadParameters(Serializable::IArchive& ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<double>(ar, "min", min_);
      loadParameter<double>(ar, "max", max_);
      loadParameter<double>(ar, "nBin", nBin_);
      ar & accumulator_;
      //loadParamComposite(in, accumulator_);

      // Validate
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (!feq(accumulator_.min(),min_)) {
         UTIL_THROW("Inconsistent values of min");
      }
      if (!feq(accumulator_.max(), max_)) {
         UTIL_THROW("Inconsistent values of max");
      }
      if (accumulator_.nBin() != nBin_) {
         UTIL_THROW("Inconsistent values of max");
      }

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void BondLengthDist::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Clear accumulator.
   */
   void BondLengthDist::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Object is not initialized"); 
      }
      accumulator_.clear();
   }
 
   /// Add particle pairs to RDF histogram.
   void BondLengthDist::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         System::MoleculeIterator molIter;
         Molecule::BondIterator   bondIter;
         double    lsq, l;

         for (system().begin(speciesId_, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(bondIter); bondIter.notEnd(); ++bondIter) {
               lsq = system().boundary().distanceSq( bondIter->atom(0).position(), 
                                                     bondIter->atom(1).position());
               l = sqrt(lsq);
               accumulator_.sample(l);
            }
         }
      }
   }  


   /// Output results to file after simulation is completed.
   void BondLengthDist::output() 
   {  

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output statistical analysis to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}
