#ifndef DDMD_CONFIG_IO_FACTORY_CPP
#define DDMD_CONFIG_IO_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigIoFactory.h"  

// Subclasses of ConfigIo 
#include "DdMdConfigIo.h"
#include "DdMdOrderedConfigIo.h"
#include "LammpsOrderedConfigIo.h"

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor
   */
   ConfigIoFactory::ConfigIoFactory(Simulation& simulation)
    : simulationPtr_(&simulation)
   {}

   /* 
   * Return a pointer to a instance of ConfigIo subclass className.
   */
   ConfigIo* ConfigIoFactory::factory(const std::string &className) const
   {
      ConfigIo *ptr = 0;

      // Try subfactories first.
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "DdMdConfigIo") {
         ptr = new DdMdConfigIo(*simulationPtr_);
      } else 
      if (className == "DdMdOrderedConfigIo") {
         ptr = new DdMdOrderedConfigIo(*simulationPtr_);
      }
      if (className == "LammpsOrderedConfigIo") {
         ptr = new LammpsOrderedConfigIo(*simulationPtr_);
      }
 
      return ptr;
   }

}
#endif
