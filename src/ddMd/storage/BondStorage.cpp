#ifndef DDMD_BOND_STORAGE_CPP
#define DDMD_BOND_STORAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BondStorage.h"
#include <ddMd/chemistry/Bond.h>

namespace DdMd
{

   /*
   * Read parameters and allocate memory.
   */
   void BondStorage::readParam(std::istream& in)
   {
      readBegin(in, "BondStorage");
      GroupStorage<2>::readParam(in);
      readEnd(in);
   }

}

#endif
