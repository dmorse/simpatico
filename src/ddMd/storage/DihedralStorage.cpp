/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DihedralStorage.h"
#include "GroupStorage.tpp"
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
