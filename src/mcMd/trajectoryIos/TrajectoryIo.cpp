#ifndef TRAJECTORY_IO_CPP
#define MCMD_TRAJECTORY_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryIo.h"
#include <mcMd/simulation/System.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor. 
   */
   TrajectoryIo::TrajectoryIo(System& system)
    : nFrames_(0),
      boundaryPtr_(&system.boundary()),
      systemPtr_(&system),
      simulationPtr_(&system.simulation())
   {}

   /* 
   * Destructor.   
   */
   TrajectoryIo::~TrajectoryIo() 
   {}

   /*
   * Read the trajectory file header.
   */
   void TrajectoryIo::readHeader(std::fstream &file)
   {
      UTIL_THROW("This TrajectoryIo class does not implement a readHeader() method.");
   }

   /*
   * Read a single frame.
   */
   void TrajectoryIo::readFrame(std::fstream& file)
   {
       UTIL_THROW("This TrajectoryIo class does not implement a readFrame() method.");
   }

   /*
   * Write the trajectory file header.
   */
   void TrajectoryIo::writeHeader(std::fstream &file)
   {
      UTIL_THROW("This TrajectoryIo class does not implement a writeHeader() method.");
   }

   /*
   * Write a single frame.
   */
   void TrajectoryIo::writeFrame(std::fstream& file)
   {
       UTIL_THROW("This TrajectoryIo class does not implement a writeFrame() method.");
   }

} 
#endif
