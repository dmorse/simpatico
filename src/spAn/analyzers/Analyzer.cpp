#ifndef SPAN_ANALYZER_CPP
#define SPAN_ANALYZER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
//#include <util/misc/FileMaster.h>
#include <util/global.h>

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   Analyzer::Analyzer(Processor& processor)
    : ParamComposite(),
      outputFileName_(),
      processorPtr_(&processor),
      interval_(1)
   {}

   /*
   * Destructor.
   */
   Analyzer::~Analyzer()
   {}

   /*
   * Read the interval from parameter file, with error checking.
   */
   void Analyzer::readInterval(std::istream &in)
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
   void Analyzer::readOutputFileName(std::istream &in)
   {  read<std::string>(in, "outputFileName", outputFileName_); }

   /*
   * Get the outputFileName string with an added suffix
   */
   std::string Analyzer::outputFileName(const std::string& suffix) const
   {
      std::string filename = outputFileName_;
      filename += suffix;
      return filename;
   }

}
#endif
