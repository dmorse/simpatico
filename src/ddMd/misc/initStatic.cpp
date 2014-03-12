#ifndef DDMD_INIT_STATIC_CPP
#define DDMD_INIT_STATIC_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "initStatic.h"
#include <ddMd/modifiers/Modifier.h>
#include <util/global.h>

namespace DdMd
{

   /*
   * Initialize all static class members in DdMd namespace.
   */
   void initStatic()
   {
      // Precondition: This function should only be called once.
      static int nCall = 0;
      if (nCall == 0) {
         // Call initStatic() methods of all relevant classes.
         Modifier::initStatic();
      }
      ++nCall;
   }

}
#endif
