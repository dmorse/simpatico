#ifndef MSDD_ATOM_STORAGE_PARAM_CPP
#define MSDD_ATOM_STORAGE_PARAM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomStorageParam.h"

namespace MsDd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   AtomStorageParam::AtomStorageParam()
    : atomCapacity_(0),
      ghostCapacity_(0),
      totalAtomCapacity_(0)
   {}
 
   /*
   * Destructor.
   */
   AtomStorageParam::~AtomStorageParam()
   {}

   /*
   * Read parameters and allocate memory.
   */
   void AtomStorageParam::readParam(std::istream& in)
   {
      readBegin(in, "AtomStorage");
      read<int>(in, "atomCapacity", atomCapacity_);
      read<int>(in, "ghostCapacity", ghostCapacity_);
      read<int>(in, "totalAtomCapacity", totalAtomCapacity_);
      readEnd(in);
   }

}
#endif
