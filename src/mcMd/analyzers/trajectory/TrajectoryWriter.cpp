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
   TrajectoryWriter::TrajectoryWriter(MdSystem& system, bool isBinary)
    : SystemAnalyzer<MdSystem>(system),
      nSample_(0),
      simulationPtr_(&system.simulation()),
      boundaryPtr_(&system.boundary()),
      isInitialized_(false),
      isBinary_(isBinary)
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
   {  ar & *this; }

   /*
   * Read interval and outputFileName.
   */
   void TrajectoryWriter::setup()
   {
      UTIL_CHECK(isInitialized_);
      std::string filename = outputFileName();
      if (isBinary()) {
         fileMaster().openOutputFile(filename, outputFile_,
                                     std::ios::out | std::ios::binary);
      } else {
         fileMaster().openOutputFile(filename, outputFile_);
      }
      nSample_ = 0;
      writeHeader();
   }

   /*
   * Dump configuration to file
   */
   void TrajectoryWriter::sample(long iStep)
   {
      UTIL_CHECK(isInitialized_);
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
