#ifndef MDPP_CONFIG_IO_CPP
#define MDPP_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigIo.h"

namespace MdPp
{

   using namespace Util;

   /*
   * Constructor.
   */
   ConfigIo::ConfigIo()
     : storagePtr_(0)
   {  setClassName("ConfigIo"); }

   /*
   * Constructor.
   */
   ConfigIo::ConfigIo(Storage& storage)
     : storagePtr_(&storage)
   {
      setClassName("ConfigIo"); 
   }

   /*
   * Destructor.
   */
   ConfigIo::~ConfigIo()
   {}

}
#endif
