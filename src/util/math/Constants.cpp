#ifndef CONSTANTS_CPP
#define CONSTANTS_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Constants.h"
#include <util/global.h>
#include <cmath>

namespace Util
{

   const double               Constants::Pi = 2.0*acos(0.0);
   const std::complex<double> Constants::Im = std::complex<double>(0.0, 1.0);

   /*
   * This function exists to guarantee initialization of static constants.
   * Calling it guarantees that the contents of this file will be linked,
   * rather than optimized away. It may only be called once.
   */
   void Constants::initStatic() 
   {
      // This method may only be called once.
      static int nCall = 0;
      ++nCall;
   }

}
#endif
