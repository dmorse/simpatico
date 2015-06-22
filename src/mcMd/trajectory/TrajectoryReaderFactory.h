#ifndef MCMD_TRAJECTORY_READER_FACTORY_H
#define MCMD_TRAJECTORY_READER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <mcMd/trajectory/TrajectoryReader.h>

#include <string>

namespace McMd
{

   using namespace Util;
   class System;

   /**
   * Default Factory for subclasses of TrajectoryReader.
   */
   class TrajectoryReaderFactory : public Factory<TrajectoryReader> 
   {

   public:

      /// Constructor
      TrajectoryReaderFactory(System& system);

      /**
      * Method to create any TrajectoryReader supplied with Simpatico.
      *
      * \param trajectoryIoName name of the TrajectoryReader subclass
      * \return TrajectoryReader* pointer to new instance of trajectryIoName
      */
      TrajectoryReader* factory(const std::string &trajectoryIoName) const;

   private:

      /// Pointer to a parent System.
      System* systemPtr_;

   };

}
#endif
