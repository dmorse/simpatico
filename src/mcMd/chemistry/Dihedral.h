#ifndef DIHEDRAL_H
#define DIHEDRAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Group.h"

namespace McMd
{

   using namespace Util;

   /**
   * Group of the 4 Atoms in a dihedral interaction.
   *
   * \ingroup Chemistry_Module
   */
   typedef Group<4> Dihedral;

}

#endif
