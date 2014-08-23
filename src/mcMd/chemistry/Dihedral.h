#ifndef MCMD_DIHEDRAL_H
#define MCMD_DIHEDRAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Group.h"

namespace McMd
{

   using namespace Util;

   /**
   * Group of the 4 Atoms in a dihedral interaction.
   *
   * \ingroup McMd_Chemistry_Module
   */
   typedef Group<4> Dihedral;

}

#endif
