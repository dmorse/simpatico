#ifndef UTIL_BIT_CPP
#define UTIL_BIT_CPP

#include "Bit.h"
#include <util/global.h>
#include <climits>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   /*
   * Constructor.
   */
   Bit::Bit(unsigned int shift)
   {  
      if (shift > sizeof(unsigned int)*CHAR_BIT) {
         UTIL_THROW("Shift is too large");
      }
      mask_ = 1 << shift;  
   }

}
#endif
