#ifndef DDMD_WRITE_TRAJECTORY_CPP
#define DDMD_WRITE_TRAJECTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "WriteTrajectory.h"
#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   WriteTrajectory::WriteTrajectory(Simulation& simulation) 
    : Analyzer(simulation),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("WriteTrajectory"); }

   /*
   * Read interval and outputFileName. 
   */
   void WriteTrajectory::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void WriteTrajectory::loadParameters(Serializable::IArchive &ar)
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
   void WriteTrajectory::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;
   }

   /*
   * Clear sample counter and close file.
   */
   void WriteTrajectory::clear() 
   {  
      nSample_ = 0; 
      if (outputFile_.is_open()) {
         outputFile_.close();
      }
   }

   /*
   * Write frame to file, header on first sample.
   */
   void WriteTrajectory::sample(long iStep) 
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
   void WriteTrajectory::output()
   {  clear(); }

}
#endif 
