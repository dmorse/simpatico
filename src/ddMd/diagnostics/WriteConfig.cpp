#ifndef DDMD_WRITE_CONFIG_CPP
#define DDMD_WRITE_CONFIG_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "WriteConfig.h"
//#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   WriteConfig::WriteConfig(Simulation& simulation) 
    : Diagnostic(simulation),
      nSample_(0),
      isInitialized_(false)
   {}

   /*
   * Read interval and outputFileName. 
   */
   void WriteConfig::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   /*
   * Read interval and outputFileName. 
   */
   void WriteConfig::clear() 
   {  nSample_ = 0; }

   /*
   * Dump configuration to file
   */
   void WriteConfig::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         // Construct new fileName: outputFileName + toString(nSample)
         std::string filename;
         filename  = outputFileName();
         filename += toString(nSample_);

         // Open output file, write data, and close file
         //simulation().fileMaster().openOutputFile(filename, outputFile_);
         simulation().writeConfig(filename);
         //outputFile_.close();
         ++nSample_;

      }
   }

}
#endif 
