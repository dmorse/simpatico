/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
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
