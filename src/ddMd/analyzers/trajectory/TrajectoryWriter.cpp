/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"
#include <ddMd/simulation/Simulation.h>
#include <util/mpi/MpiLoader.h>

#include <fstream>
#include <ios>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   TrajectoryWriter::TrajectoryWriter(Simulation& simulation, bool isBinary)
    : Analyzer(simulation),
      isInitialized_(false),
      isBinary_(isBinary),
      domainPtr_(0),
      boundaryPtr_(0),
      atomStoragePtr_(0)
      #ifdef SIMP_BOND
      , bondStoragePtr_(0)
      #endif
      #ifdef SIMP_ANGLE
      , angleStoragePtr_(0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , dihedralStoragePtr_(0)
      #endif
   {
      setClassName("TrajectoryWriter");
      domainPtr_ = &simulation.domain();
      boundaryPtr_ = &simulation.boundary();

      atomStoragePtr_ = &simulation.atomStorage();
      #ifdef SIMP_BOND
      bondStoragePtr_ = &simulation.bondStorage();
      #endif
      #ifdef SIMP_ANGLE
      angleStoragePtr_ = &simulation.angleStorage();
      #endif
      #ifdef SIMP_DIHEDRAL
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
      isInitialized_ = true;
   }

   /*
   * Save internal state to output archive.
   */
   void TrajectoryWriter::save(Serializable::OArchive& ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
   }

   /*
   * Setup - open the trajectory file.
   */
   void TrajectoryWriter::setup()
   {  
      FileMaster& fileMaster = simulation().fileMaster();
      if (isIoProcessor()) {
         if (isBinary()) {
            fileMaster.openOutputFile(outputFileName(), outputFile_, 
                                      std::ios::out | std::ios::binary);
         } else {
            fileMaster.openOutputFile(outputFileName(), outputFile_);
         }
      }
      writeHeader(outputFile_);
   }

   /*
   * Write a frame to file.
   */
   void TrajectoryWriter::sample(long iStep)
   {
      if (isAtInterval(iStep))  {
         writeFrame(outputFile_, iStep);
      }
   }

   /*
   * Clear sample counter and close file.
   */
   void TrajectoryWriter::clear()
   {
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
