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
      fileMasterPtr_(0),
      interval_(1)
   {}
   
   /* 
   * Destructor.
   */
   Diagnostic::~Diagnostic()
   {}
   
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

   /*
   * Load interval from archive.
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
   
      // Read interval value (inherited from Interval)
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
   * Save interval to an archive.
   */
   void Diagnostic::saveInterval(Serializable::OArchive& ar)
   {  ar & interval_; }

   /*
   * Read output file name and open output file.
   */
   void Diagnostic::readOutputFileName(std::istream &in) 
   {  read<std::string>(in, "outputFileName", outputFileName_); }

   /*
   * Load output file name from archive.
   */
   void Diagnostic::loadOutputFileName(Serializable::IArchive& ar)
   {  ar & outputFileName_; }

   /*
   * Save output file name to an archive.
   */
   void Diagnostic::saveOutputFileName(Serializable::OArchive& ar)
   {  ar & outputFileName_; }

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
