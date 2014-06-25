#ifndef DDMD_SP_ANALYZER_CPP
#define DDMD_SP_ANALYZER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SpAnalyzer.h"
//#include <util/misc/FileMaster.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   SpAnalyzer::SpAnalyzer(Processor& processor)
    : ParamComposite(),
      outputFileName_(),
      processorPtr_(&processor),
      interval_(1)
   {}

   /*
   * Destructor.
   */
   SpAnalyzer::~SpAnalyzer()
   {}

   /*
   * Read the interval from parameter file, with error checking.
   */
   void SpAnalyzer::readInterval(std::istream &in)
   {

      // Read interval value (inherited from Interval)
      read<long>(in, "interval", interval_);

      // Check that interval has a nonzero, positive value
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }

   }

   /*
   * Read output file name and open output file.
   */
   void SpAnalyzer::readOutputFileName(std::istream &in)
   {  read<std::string>(in, "outputFileName", outputFileName_); }

   /*
   * Get the outputFileName string with an added suffix
   */
   std::string SpAnalyzer::outputFileName(const std::string& suffix) const
   {
      std::string filename = outputFileName_;
      filename += suffix;
      return filename;
   }

}
#endif
