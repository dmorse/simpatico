#ifndef MCMD_DUMP_CONFIG_CPP
#define MCMD_DUMP_CONFIG_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DumpConfig.h"
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DumpConfig::DumpConfig(System& system) 
    : SystemDiagnostic<System>(system),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("DumpConfig"); }

   /*
   * Read interval and outputFileName. 
   */
   void DumpConfig::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void DumpConfig::loadParameters(Serializable::IArchive& ar)
   {
      Diagnostic::loadParameters(ar);
      ar & nSample_;
      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void DumpConfig::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Read interval and outputFileName. 
   */
   void DumpConfig::setup() 
   {  nSample_ = 0; }

   /*
   * Dump configuration to file
   */
   void DumpConfig::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         // Construct new fileName: outputFileName + toString(nSample)
         std::string filename;
         filename  = outputFileName();
         filename += toString(nSample_);

         // Open output file, write data, and close file
         fileMaster().openOutputFile(filename, outputFile_);
         system().writeConfig(outputFile_);
         outputFile_.close();
         ++nSample_;

      }
   }
  
}
#endif // DUMP_CONFIG_CPP
