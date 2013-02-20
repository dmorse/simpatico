#ifndef MCMD_TRAJECTORY_IO_FACTORY_CPP
#define MCMD_TRAJECTORY_IO_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryIoFactory.h"

// Subclasses of ConfigIo
#include "DCDTrajectoryIo.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor
   */
   TrajectoryIoFactory::TrajectoryIoFactory(System& system)
    : systemPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of Trajectory subclass className.
   */
   TrajectoryIo* TrajectoryIoFactory::factory(const std::string &className) const
   {
      TrajectoryIo *ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "DCDTrajectoryIo") {
         ptr = new DCDTrajectoryIo(*systemPtr_);
      }
      return ptr;
   }

}
#endif
