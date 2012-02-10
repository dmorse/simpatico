#ifndef PLAN_CPP
#define PLAN_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Plan.h"

namespace DdMd
{

   unsigned int Plan::GMask[3][2] = { {0x0001, 0x0002}, {0x0004, 0x0008}, {0x0010, 0x0020} };
   unsigned int Plan::EMask[3][2] = { {0x0100, 0x0200}, {0x0400, 0x0800}, {0x1000, 0x2000} };

}
#endif
