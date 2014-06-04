#ifndef DDMD_SP_CONFIG_IO_CPP
#define DDMD_SP_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SpConfigIo.h"

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   SpConfigIo::SpConfigIo()
     : configurationPtr_(0)
   {  setClassName("SpConfigIo"); }

   /*
   * Constructor.
   */
   SpConfigIo::SpConfigIo(SpConfiguration& configuration)
     : configurationPtr_(&configuration)
   {
      setClassName("SpConfigIo"); 
   }

   /*
   * Destructor.
   */
   SpConfigIo::~SpConfigIo()
   {}

}
#endif
