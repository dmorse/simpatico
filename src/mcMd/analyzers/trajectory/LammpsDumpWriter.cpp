/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsDumpWriter.h"
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
   LammpsDumpWriter::LammpsDumpWriter(System& system) 
    : SystemAnalyzer<System>(system),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("LammpsDumpWriter"); }

   /*
   * Read interval and outputFileName. 
   */
   void LammpsDumpWriter::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void LammpsDumpWriter::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);
      ar & nSample_;
      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void LammpsDumpWriter::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Read interval and outputFileName. 
   */
   void LammpsDumpWriter::setup() 
   {  
      nSample_ = 0; 
      std::string filename;
      filename  = outputFileName();
      fileMaster().openOutputFile(filename, outputFile_);
   }

   /*
   * Dump configuration to file
   */
   void LammpsDumpWriter::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         system().writeConfig(outputFile_);
         ++nSample_;
      }
   }
  
   /*
   * Close output file. 
   */
   void LammpsDumpWriter::output() 
   {  outputFile_.close(); }

}
