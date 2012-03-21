#ifndef MCMD_TRAJECTORY_IO_FACTORY_H
#define MCMD_TRAJECTORY_IO_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <mcMd/trajectoryIos/TrajectoryIo.h>

#include <string>

namespace McMd
{

   using namespace Util;
   class System;

   /**
   * Default Factory for subclasses of TrajectoryIo.
   */
   class TrajectoryIoFactory : public Factory<TrajectoryIo> 
   {

   public:

      /// Constructor
      TrajectoryIoFactory(System& system);

      /**
      * Method to create any TrajectoryIo supplied with Simpatico.
      *
      * \param trajectoryIoName name of the TrajectoryIo subclass
      * \return TrajectoryIo* pointer to new instance of trajectryIoName
      */
      TrajectoryIo* factory(const std::string &trajectoryIoName) const;

   private:

      /// Pointer to a parent System.
      System* systemPtr_;

   };

}
#endif
