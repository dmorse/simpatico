#ifndef MSDD_BOND_STORAGE_PARAM_CPP
#define MSDD_BOND_STORAGE_PARAM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BondStorageParam.h"

namespace MsDd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   BondStorageParam::BondStorageParam()
   {}
 
   /*
   * Destructor.
   */
   BondStorageParam::~BondStorageParam()
   {}

   /*
   * Read parameters and allocate memory.
   */
   void BondStorageParam::readParam(std::istream& in)
   {
      readBegin(in, "BondStorage");
      GroupStorageParam<2>::readParam(in);
      readEnd(in);
   }

}
#endif
