#ifndef DDMD_INIT_STATIC_H
#define DDMD_INIT_STATIC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   /**
   * Guarantee initialization of all static class members in DdMd namespace.
   */
   void initStatic();

}
#endif
