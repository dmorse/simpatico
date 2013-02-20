#ifndef DDMD_DIHEDRAL_STORAGE_CPP
#define DDMD_DIHEDRAL_STORAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DihedralStorage.h"
#include <ddMd/chemistry/Dihedral.h>

namespace DdMd
{

   /*
   * Constructor  
   */
   DihedralStorage::DihedralStorage()
   { setClassName("DihedralStorage"); }

   /*
   * Read parameters and allocate memory.
   */
   void DihedralStorage::readParameters(std::istream& in)
   {
      GroupStorage<4>::readParameters(in);
   }

}
#endif
