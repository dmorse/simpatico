#ifndef DDMD_MODIFIER_CPP
#define DDMD_MODIFIER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Modifier.h"

namespace DdMd
{

   using namespace Util;

   // Constructor.
   Modifier::Modifier(Simulation& simulation)
    : simulationPtr_(&simulation),
      interval_(1)
   {}

   // Destructor.
   virtual ~Modifier()
   {};

}
#endif
