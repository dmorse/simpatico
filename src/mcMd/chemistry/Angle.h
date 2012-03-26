#ifndef MCMD_ANGLE_H
#define MCMD_ANGLE_H

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
   * Group of 3 Atoms involved in an angle potential.
   *
   * \ingroup McMd_Chemistry_Module
   */
   typedef Group<3> Angle;

}

#endif
