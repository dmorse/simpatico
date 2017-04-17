/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <tools/processor/Processor.h>
#include <util/misc/FileMaster.h>
#include <util/global.h>

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   Analyzer::Analyzer(Configuration& configuration)
    : ParamComposite(),
      outputFileName_(),
      configurationPtr_(&configuration),
      fileMasterPtr_(0),
      interval_(1),
      ownsFileMaster_(false)
   {}

   /*
   * Constructor.
   */
   Analyzer::Analyzer(Processor& processor)
    : ParamComposite(),
      outputFileName_(),
      configurationPtr_(&processor),
      fileMasterPtr_(&processor.fileMaster()),
      interval_(1),
      ownsFileMaster_(false)
   {}

   /*
   * Constructor.
   */
   Analyzer::Analyzer(Configuration& configuration, FileMaster& fileMaster)
    : ParamComposite(),
      outputFileName_(),
      configurationPtr_(&configuration),
      fileMasterPtr_(&fileMaster),
      interval_(1),
      ownsFileMaster_(false)
   {}

   /*
   * Destructor.
   */
   Analyzer::~Analyzer()
   {
      if (fileMasterPtr_) {
         if (ownsFileMaster_) {
            delete fileMasterPtr_;
         }
      }
   }

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

   /*
   * Get an associated FileMaster by reference.
   */
   FileMaster& Analyzer::fileMaster()
   {
      if (fileMasterPtr_ == 0) {
         fileMasterPtr_ = new FileMaster();
         ownsFileMaster_ = true;
      }
      return *fileMasterPtr_; 
   }

}
