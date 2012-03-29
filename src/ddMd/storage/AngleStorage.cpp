#ifndef DDMD_ANGLE_STORAGE_CPP
#define DDMD_ANGLE_STORAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AngleStorage.h"
#include <ddMd/chemistry/Angle.h>

namespace DdMd
{

   /*
   * Read parameters and allocate memory.
   */
   void AngleStorage::readParam(std::istream& in)
   {
      readBegin(in, "AngleStorage");
      GroupStorage<3>::readParam(in);
      readEnd(in);
   }

}

#endif
