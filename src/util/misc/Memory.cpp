#ifndef UTIL_MEMORY_CPP
#define UTIL_MEMORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Memory.h"

namespace Util
{

   /// Total amount of memory allocated, in bytes.
   size_t Memory::total_ = 0;

   /**
   * Call this to ensure compilation of this file. 
   */
   void Memory::initStatic()
   {   Memory::total_ = 0; }  

   /**
   * Return total amount of memory allocated thus far.
   */
   size_t Memory::total()
   { return total_; }

} 
#endif
