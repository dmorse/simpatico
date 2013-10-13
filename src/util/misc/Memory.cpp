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
   int Memory::total_ = 0;

   /// Maximum of total over course of simulation.
   int Memory::max_ = 0;

   /*
   * Call this to ensure compilation of this file. 
   */
   void Memory::initStatic()
   {   
      Memory::total_ = 0; 
      Memory::max_ = 0; 
   }  

   /*
   * Return total amount of memory allocated thus far.
   */
   int Memory::total()
   { return total_; }

   /*
   * Return maximum amount of allocated memory thus far.
   */
   int Memory::max()
   { return max_; }

   /*
   * Compute memory usage statistics (call on all processors).
   */
   #ifdef UTIL_MPI
   int Memory::max(MPI::Intracomm& communicator)
   { 
      int maxGlobal;
      int maxLocal = max_;
      communicator.Allreduce(&maxLocal, &maxGlobal, 1, MPI::INT, MPI::MAX);
      return maxGlobal;
   }
   #endif

} 
#endif
