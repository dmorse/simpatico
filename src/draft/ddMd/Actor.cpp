#ifndef DDMD_ACTOR_CPP
#define DDMD_ACTOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Actor.h"

namespace DdMd
{

   using namespace Util;

   // Constructor.
   Actor::Actor(Simulation& simulation)
    : interval_(1),
      hasSetupPostExchange_(false),
      hasSetupPostNeighbor_(false),
      hasSetupPostForce_(false),
      hasPreIntegrate_(false),
      hasPostIntegrate_(false),
      hasPreTransform_(false),
      hasPreExchange_(false),
      hasPostExchange_(false),
      hasPostNeighbor_(false),
      hasPreUpdate_(false),
      hasPostUpdate_(false),
      hasPreForce_(false),
      hasPostForce_(false),
      hasEndOfStep_(false),
      hasPackExchange_(false),
      hasUnpackExchange_(false),
      hasPackUpdate_(false),
      hasUnpackUpdate_(false),
      hasPackReverseUpdate_(false),
      hasUnpackReverseUpdate_(false),
      simulationPtr_(&simulation)
   {}

   // Destructor.
   Actor::~Actor()
   {}

}
#endif
