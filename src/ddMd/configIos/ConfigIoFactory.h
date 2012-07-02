#ifndef DDMD_CONFIG_IO_FACTORY_H
#define DDMD_CONFIG_IO_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <ddMd/configIos/ConfigIo.h>

#include <string>

namespace DdMd
{

   using namespace Util;
   class Simulation;

   /**
   * Default Factory for subclasses of ConfigIo.
   */
   class ConfigIoFactory : public Factory<ConfigIo> 
   {

   public:

      /// Constructor
      ConfigIoFactory(Simulation& simulation);

      /**
      * Method to create any species supplied with Simpatico.
      *
      * \param  subclassName  name of the ConfigIo subclass
      * \return ConfigIo*     pointer to new instance of subclassName
      */
      ConfigIo* factory(const std::string &subclassName) const;

   private:

      /// Pointer to a parent Simulation.
      Simulation* simulationPtr_;

   };

}
#endif
