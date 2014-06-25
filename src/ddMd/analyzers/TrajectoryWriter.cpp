#ifndef DDMD_TRAJECTORY_WRITER_CPP
#define DDMD_TRAJECTORY_WRITER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"
#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   TrajectoryWriter::TrajectoryWriter(Simulation& simulation) 
    : Analyzer(simulation),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("TrajectoryWriter"); }

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
   * Load internal state from an archive.
   */
   void TrajectoryWriter::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void TrajectoryWriter::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;
   }

   /*
   * Clear sample counter and close file.
   */
   void TrajectoryWriter::clear() 
   {  
      nSample_ = 0; 
      if (outputFile_.is_open()) {
         outputFile_.close();
      }
   }

   /*
   * Write frame to file, header on first sample.
   */
   void TrajectoryWriter::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         if (nSample_ == 0) {
            simulation().fileMaster().openOutputFile(outputFileName(), outputFile_);
            writeHeader(outputFile_, iStep);
         }
         writeFrame(outputFile_, iStep);
         ++nSample_;
      }
   }

   /*
   * Clear sample counter and close output file.
   */
   void TrajectoryWriter::output()
   {  clear(); }

}
#endif 
