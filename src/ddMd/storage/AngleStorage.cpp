#ifndef DDMD_ANGLE_STORAGE_CPP
#define DDMD_ANGLE_STORAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AngleStorage.h"
#include "GroupStorage.tpp"
#include <ddMd/chemistry/Angle.h>

namespace DdMd
{

   /*
   * Constructor  
   */
   AngleStorage::AngleStorage()
   { setClassName("AngleStorage"); }

   /*
   * Read parameters and allocate memory.
   */
   void AngleStorage::readParameters(std::istream& in)
   {
      GroupStorage<3>::readParameters(in);
   }

}

#endif
