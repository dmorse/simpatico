/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigWriter.h"
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
   ConfigWriter::ConfigWriter(System& system) 
    : SystemAnalyzer<System>(system),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("ConfigWriter"); }

   /*
   * Read interval and outputFileName. 
   */
   void ConfigWriter::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void ConfigWriter::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);
      ar & nSample_;
      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void ConfigWriter::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Read interval and outputFileName. 
   */
   void ConfigWriter::setup() 
   {  nSample_ = 0; }

   /*
   * Dump configuration to file
   */
   void ConfigWriter::sample(long iStep) 
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
