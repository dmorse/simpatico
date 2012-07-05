#ifndef MCMD_LOG_PROGRESS_CPP
#define MCMD_LOG_PROGRESS_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "LogProgress.h"
#include <util/util/Log.h>

namespace McMd
{

   using namespace Util;

   /*
   * Read interval.
   */
   void LogProgress::readParam(std::istream& in) 
   {  readInterval(in); }

   /*
   * Dump configuration to file
   */
   void LogProgress::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Log::file() << "iStep = " << iStep << std::endl;
      }
   }
  
}
#endif 
