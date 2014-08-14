#ifndef SPAN_CONFIG_IO_CPP
#define SPAN_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigIo.h"

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   ConfigIo::ConfigIo()
     : configurationPtr_(0)
   {  setClassName("ConfigIo"); }

   /*
   * Constructor.
   */
   ConfigIo::ConfigIo(Configuration& configuration)
     : configurationPtr_(&configuration)
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
