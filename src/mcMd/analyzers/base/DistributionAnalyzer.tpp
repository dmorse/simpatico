/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DistributionAnalyzer.h"
#include <util/math/feq.h>
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   template <class SystemType>
   DistributionAnalyzer<SystemType>::DistributionAnalyzer(SystemType& system) 
    : SystemAnalyzer<SystemType>(system),
      min_(0.0),
      max_(0.0),
      nBin_(0)
   {  setClassName("DistributionAnalyzer"); }

   /*
   * Read parameters from file, set up and clear accumulator
   */
   template <class SystemType>
   void DistributionAnalyzer<SystemType>::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      ParamComposite::template read<double>(in, "min", min_);
      ParamComposite::template read<double>(in, "max", max_);
      ParamComposite::template read<double>(in, "nBin", nBin_);
      accumulator_.setParam(min_, max_, nBin_);
      accumulator_.clear();
   }

   /*
   * Load state from an archive.
   */
   template <class SystemType>
   void 
   DistributionAnalyzer<SystemType>::loadParameters(Serializable::IArchive& ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      ParamComposite::template loadParameter<double>(ar, "min", min_);
      ParamComposite::template loadParameter<double>(ar, "max", max_);
      ParamComposite::template loadParameter<double>(ar, "nBin", nBin_);
      ar & accumulator_;

      // Validate
      if (!feq(accumulator_.min(),min_)) {
         UTIL_THROW("Inconsistent values of min");
      }
      if (!feq(accumulator_.max(), max_)) {
         UTIL_THROW("Inconsistent values of max");
      }
      if (accumulator_.nBin() != nBin_) {
         UTIL_THROW("Inconsistent values of max");
      }

   }

   /*
   * Save state to archive.
   */
   template <class SystemType>
   void DistributionAnalyzer<SystemType>::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Clear accumulator.
   */
   template <class SystemType>
   void DistributionAnalyzer<SystemType>::setup() 
   {
      accumulator_.clear();
   }
 
   /*
   * Add value to histogram.
   */
   template <class SystemType>
   void DistributionAnalyzer<SystemType>::increment(double value) 
   {
      accumulator_.sample(value);
   }  

   /*
   * Output results to file after simulation is completed.
   */
   template <class SystemType>
   void DistributionAnalyzer<SystemType>::output() 
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
