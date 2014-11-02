#ifndef DDMD_ANALYZER_CPP
#define DDMD_ANALYZER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <util/misc/FileMaster.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   // Define and initialize static variable
   long Analyzer::baseInterval = 1;

   void Analyzer::initStatic()
   {  Analyzer::baseInterval = 1; }

   /*
   * Default constructor.
   */
   Analyzer::Analyzer(Simulation& simulation)
    : ParamComposite(),
      outputFileName_(),
      simulationPtr_(&simulation),
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

      // Check that baseInterval has a nonzero, positive value
      if (baseInterval == 0) {
         UTIL_THROW("baseInterval == 0");
      }
      if (baseInterval < 0) {
         UTIL_THROW("baseInterval < 0");
      }

      // Read interval value (inherited from Interval)
      read<long>(in, "interval", interval_);

      // Check that interval has a nonzero, positive value
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }

      // Check that interval is a multiple of baseInterval
      if (interval_ % baseInterval != 0) {
         UTIL_THROW("interval is not a multiple of baseInterval");
      }

   }

   /*
   * Load parameter interval from input archive.
   */
   void Analyzer::loadInterval(Serializable::IArchive &ar)
   {
      // Check that baseInterval has a nonzero, positive value
      if (baseInterval == 0) {
         UTIL_THROW("baseInterval == 0");
      }
      if (baseInterval < 0) {
         UTIL_THROW("baseInterval < 0");
      }

      loadParameter<long>(ar, "interval", interval_);

      // Check that interval has a nonzero, positive value
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }

      // Check that interval is a multiple of baseInterval
      if (interval_ % baseInterval != 0) {
         UTIL_THROW("interval is not a multiple of baseInterval");
      }
   }

   /*
   * Save interval parameter to an archive.
   */
   void Analyzer::saveInterval(Serializable::OArchive &ar)
   {  ar << interval_; }

   /*
   * Read output file name and open output file.
   */
   void Analyzer::readOutputFileName(std::istream &in)
   {  read<std::string>(in, "outputFileName", outputFileName_); }

   /*
   * Load output file name to an archive.
   */
   void Analyzer::loadOutputFileName(Serializable::IArchive &ar)
   {
      loadParameter<std::string>(ar, "outputFileName", outputFileName_);
   }

   /*
   * Save output file name to an archive.
   */
   void Analyzer::saveOutputFileName(Serializable::OArchive &ar)
   {  ar << outputFileName_; }

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
