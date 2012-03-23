#ifndef MEMORY_CPP
#define MEMORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Memory.h"

namespace Util
{

   size_t Memory::totalAllocated_ = 0;

   /**
   * Call this to ensure compilation of this file. 
   */
   void Memory::initStatic()
   {   Memory::totalAllocated_ = 0; }  

   /**
   * Return total amount of memory allocated thus far.
   */
   size_t Memory::totalAllocated()
   { return totalAllocated_; }

} 
#endif
