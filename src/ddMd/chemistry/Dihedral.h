#ifndef DDMD_DIHEDRAL_H
#define DDMD_DIHEDRAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Group.h"

namespace DdMd
{

   typedef Group<4> Dihedral;

   #if 0
   /**
   * A group of four atoms interacting via an angle potential.
   *
   * \ingroup DdMd_Chemistry_Module
   */
   class Dihedral : public Group<4>
   {};
   #endif

}
#endif
