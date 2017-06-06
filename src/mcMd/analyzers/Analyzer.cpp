/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   // Define and initialize static variable
   long Analyzer::baseInterval = 0;

   void Analyzer::initStatic()
   {  Analyzer::baseInterval = 0; }

   /* 
   * Default constructor.
   */
   Analyzer::Analyzer()
    : ParamComposite(),
      outputFileName_(),
      interval_(1),
      fileMasterPtr_(0)
   {}
   
   /* 
   * Destructor.
   */
   Analyzer::~Analyzer()
   {}
   
   /*
   * Read parameters from stream, default implementation.
   */
   void Analyzer::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
   }

   /*
   * Read the interval from parameter file, with error checking.
   */
   void Analyzer::readInterval(std::istream &in, bool checkBase) 
   {
      if (checkBase) {
         // Require that baseInterval is positive
         if (baseInterval == 0) {
            UTIL_THROW("baseInterval == 0");
         }
         if (baseInterval < 0) {
            UTIL_THROW("baseInterval < 0");
         }
      }
   
      // Read interval value (inherited from Interval)
      read<long>(in, "interval", interval_);
   
      // Postconditons
      // Require that interval is positive
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }
      if (checkBase) {
         if (interval_ % baseInterval != 0) {
            UTIL_THROW("interval not multiple of baseInterval");
         }
      }

   }

   void Analyzer::readOutputFileName(std::istream &in) 
   {  read<std::string>(in, "outputFileName", outputFileName_); }

   /*
   * Load parameters from archive, default implementation.
   */
   void Analyzer::loadParameters(Serializable::IArchive& ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
   }

   /*
   * Load parameters from archive, with error checking.
   */
   void Analyzer::loadInterval(Serializable::IArchive& ar,
                               bool checkBase)
   {
      if (checkBase) {
         // Check that Analyzer::baseInterval is positive
         if (baseInterval == 0) {
            UTIL_THROW("baseInterval == 0");
         }
         if (baseInterval < 0) {
            UTIL_THROW("baseInterval < 0");
         }
      }
   
      // Load parameters
      loadParameter<long>(ar, "interval", interval_);
   
      // Postconditons
      // Require that interval is positive integer 
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }
      if (checkBase) {
         if (interval_ % baseInterval != 0) {
            UTIL_THROW("interval not multiple of baseInterval");
         }
      }

   }

   /*
   * Load outputFileName from archive.
   */
   void Analyzer::loadOutputFileName(Serializable::IArchive& ar)
   { loadParameter<std::string>(ar, "outputFileName", outputFileName_); }

   /*
   * Save interval and outputFileName to archive.
   */
   void Analyzer::save(Serializable::OArchive& ar)
   {
      ar & interval_;
      ar & outputFileName_;
   }

   /*
   * Set the FileMaster.
   */
   void Analyzer::setFileMaster(FileMaster& fileMaster)
   {  fileMasterPtr_ = &fileMaster; }

   /*
   * Get the FileMaster by reference.
   */
   FileMaster& Analyzer::fileMaster()
   {  
      assert(fileMasterPtr_);
      return (*fileMasterPtr_);
   }

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
