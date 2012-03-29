#ifndef DDMD_DIHEDRAL_STORAGE_CPP
#define DDMD_DIHEDRAL_STORAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DihedralStorage.h"
#include <ddMd/chemistry/Dihedral.h>

namespace DdMd
{

   /*
   * Read parameters and allocate memory.
   */
   void DihedralStorage::readParam(std::istream& in)
   {
      readBegin(in, "DihedralStorage");
      GroupStorage<4>::readParam(in);
      readEnd(in);
   }

}
#endif
