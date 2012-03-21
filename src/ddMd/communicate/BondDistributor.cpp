#ifndef DDMD_BOND_DISTRIBUTOR_CPP
#define DDMD_BOND_DISTRIBUTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BondDistributor.h"
#include "GroupDistributor.cpp"

namespace DdMd
{

   /*
   * Read cacheCapacity and allocate all required memory.
   */
   void BondDistributor::readParam(std::istream& in)
   {
      // Read parameter file block
      readBegin(in, "BondDistributor");
      read<int>(in, "cacheCapacity", cacheCapacity_);
      readEnd(in);
 
      // Do actual allocation
      allocate();
   }

}
#endif
