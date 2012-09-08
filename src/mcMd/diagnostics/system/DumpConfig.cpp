#ifndef MCMD_DUMP_CONFIG_CPP
#define MCMD_DUMP_CONFIG_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DumpConfig.h"
#include <mcMd/util/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/util/ioUtil.h>

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
  
   /*
   * Save state to binary file archive.
   */
   void DumpConfig::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void DumpConfig::load(Serializable::IArchiveType& ar)
   { ar & *this; }

}
#endif // DUMP_CONFIG_CPP
