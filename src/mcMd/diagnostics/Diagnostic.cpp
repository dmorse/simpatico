#ifndef MCMD_DIAGNOSTIC_CPP
#define MCMD_DIAGNOSTIC_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Diagnostic.h"
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   // Define and initialize static variable
   long Diagnostic::baseInterval = 1;

   void Diagnostic::initStatic()
   {  Diagnostic::baseInterval = 1; }

   /* 
   * Default constructor.
   */
   Diagnostic::Diagnostic()
    : ParamComposite(),
      outputFileName_(),
      interval_(1),
      fileMasterPtr_(0)
   {}
   
   /* 
   * Destructor.
   */
   Diagnostic::~Diagnostic()
   {}
   
   /*
   * Read parameters from stream, default implementation.
   */
   void Diagnostic::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
   }

   /*
   * Read the interval from parameter file, with error checking.
   */
   void Diagnostic::readInterval(std::istream &in) 
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
   
      // Postconditons
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }
      if (interval_ % baseInterval != 0) {
         UTIL_THROW("interval is not a multiple of baseInterval");
      }
   }

   void Diagnostic::readOutputFileName(std::istream &in) 
   {  read<std::string>(in, "outputFileName", outputFileName_); }

   /*
   * Load parameters from archive, default implementation.
   */
   void Diagnostic::loadParameters(Serializable::IArchive& ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
   }

   /*
   * Load parameters from archive, with error checking.
   */
   void Diagnostic::loadInterval(Serializable::IArchive& ar)
   {
      // Check that baseInterval has a nonzero, positive value
      if (baseInterval == 0) {
         UTIL_THROW("baseInterval == 0");
      }
      if (baseInterval < 0) {
         UTIL_THROW("baseInterval < 0");
      }
   
      // Load parameters
      loadParameter<long>(ar, "interval", interval_);
   
      // Postconditons
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }
      if (interval_ % baseInterval != 0) {
         UTIL_THROW("interval is not a multiple of baseInterval");
      }
   }

   /*
   * Load outputFileName from archive.
   */
   void Diagnostic::loadOutputFileName(Serializable::IArchive& ar)
   { loadParameter<std::string>(ar, "outputFileName", outputFileName_); }

   /*
   * Save interval and outputFileName to archive.
   */
   void Diagnostic::save(Serializable::OArchive& ar)
   {
      ar & interval_;
      ar & outputFileName_;
   }

   /*
   * Set the FileMaster.
   */
   void Diagnostic::setFileMaster(FileMaster& fileMaster)
   {  fileMasterPtr_ = &fileMaster; }

   /*
   * Get the FileMaster by reference.
   */
   FileMaster& Diagnostic::fileMaster()
   {  
      assert(fileMasterPtr_);
      return (*fileMasterPtr_);
   }

   /*
   * Get the outputFileName string with an added suffix
   */
   std::string Diagnostic::outputFileName(const std::string& suffix) const
   {
      std::string filename = outputFileName_;
      filename += suffix;
      return filename;
   }

}
#endif
