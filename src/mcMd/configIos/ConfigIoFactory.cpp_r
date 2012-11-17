#ifndef MCMD_CONFIG_IO_FACTORY_CPP
#define MCMD_CONFIG_IO_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigIoFactory.h"  

// Subclasses of ConfigIo 
#include "McConfigIo.h"
#include "MdConfigIo.h"
#include "DdMdConfigIo.h"
#include "LammpsConfigIo.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor
   */
   ConfigIoFactory::ConfigIoFactory(System& system)
    : systemPtr_(&system)
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

      if (className == "McConfigIo") {
         ptr = new McConfigIo(*systemPtr_);
      } else
      if (className == "MdConfigIo") {
         ptr = new MdConfigIo(*systemPtr_);
      } else 
      if (className == "LammpsConfigIo") {
         ptr = new LammpsConfigIo(*systemPtr_);
      } else
      if (className == "DdMdConfigIo") {
         ptr = new DdMdConfigIo(*systemPtr_);
      }
      return ptr;
   }

}
#endif
