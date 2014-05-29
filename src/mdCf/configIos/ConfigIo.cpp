#ifndef MDCF_CONFIG_IO_CPP
#define MDCF_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigIo.h"

namespace MdCf
{

   using namespace Util;

   /*
   * Constructor.
   */
   ConfigIo::ConfigIo()
     : systemPtr_(0)
   {  setClassName("ConfigIo"); }

   /*
   * Constructor.
   */
   ConfigIo::ConfigIo(System& system)
     : systemPtr_(&system)
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
