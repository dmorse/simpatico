#ifndef CONFIG_IO_FACTORY_H
#define CONFIG_IO_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <mcMd/configIos/ConfigIo.h>

#include <string>

namespace McMd
{

   using namespace Util;
   class System;

   /**
   * Default Factory for subclasses of ConfigIo.
   */
   class ConfigIoFactory : public Factory<ConfigIo> 
   {

   public:

      /// Constructor
      ConfigIoFactory(System& system);

      /**
      * Method to create any species supplied with Simpatico.
      *
      * \param speciesName name of the ConfigIo subclass
      * \return ConfigIo* pointer to new instance of speciesName
      */
      ConfigIo* factory(const std::string &speciesName) const;

   private:

      /// Pointer to a parent System.
      System* systemPtr_;

   };

}
#endif
