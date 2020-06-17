/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"
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
   TrajectoryWriter::TrajectoryWriter(System& system) 
    : SystemAnalyzer<System>(system),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("TrajectoryWriter"); }

  
   /*
   * Destructor.
   */
   TrajectoryWriter::~TrajectoryWriter() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   void TrajectoryWriter::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void TrajectoryWriter::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);
      ar & nSample_;
      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void TrajectoryWriter::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Read interval and outputFileName. 
   */
   void TrajectoryWriter::setup() 
   {  
      nSample_ = 0; 
      std::string filename;
      filename  = outputFileName();
      fileMaster().openOutputFile(filename, outputFile_);
      writeHeader();
   }

   /*
   * Dump configuration to file
   */
   void TrajectoryWriter::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {
         writeFrame(iStep);
         ++nSample_;
      }
   }
  
   /*
   * Read interval and outputFileName. 
   */
   void TrajectoryWriter::output() 
   {  outputFile_.close(); }

}
