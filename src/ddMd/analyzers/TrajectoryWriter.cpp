#ifndef DDMD_TRAJECTORY_WRITER_CPP
#define DDMD_TRAJECTORY_WRITER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"
#include <ddMd/simulation/Simulation.h>
#include <util/mpi/MpiLoader.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   TrajectoryWriter::TrajectoryWriter(Simulation& simulation)
    : Analyzer(simulation),
      nSample_(0),
      isInitialized_(false),
      domainPtr_(0),
      boundaryPtr_(0),
      atomStoragePtr_(0)
      #ifdef INTER_BOND
      , bondStoragePtr_(0)
      #endif
      #ifdef INTER_ANGLE
      , angleStoragePtr_(0)
      #endif
      #ifdef INTER_DIHEDRAL
      , dihedralStoragePtr_(0)
      #endif
   {
      setClassName("TrajectoryWriter");
      domainPtr_ = &simulation.domain();
      boundaryPtr_ = &simulation.boundary();

      atomStoragePtr_ = &simulation.atomStorage();
      #ifdef INTER_BOND
      bondStoragePtr_ = &simulation.bondStorage();
      #endif
      #ifdef INTER_ANGLE
      angleStoragePtr_ = &simulation.angleStorage();
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_ = &simulation.dihedralStorage();
      #endif
   }

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
   * Save internal state to output archive.
   */
   void TrajectoryWriter::save(Serializable::OArchive& ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;
   }

   /*
   * Setup - open the trajectory file.
   */
   void TrajectoryWriter::setup()
   {  simulation().fileMaster().openOutputFile(outputFileName(), outputFile_); }

   /*
   * Write frame to file, header on first sample.
   */
   void TrajectoryWriter::sample(long iStep)
   {
      if (isAtInterval(iStep))  {
         if (nSample_ == 0) {
            writeHeader(outputFile_, iStep);
         }
         writeFrame(outputFile_, iStep);
         ++nSample_;
      }
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
   * Clear sample counter and close output file.
   */
   void TrajectoryWriter::output()
   {  clear(); }

}
#endif
