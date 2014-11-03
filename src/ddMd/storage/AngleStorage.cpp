/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
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
