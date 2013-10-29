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

   const Bit Actor::Flags::SetupExchange = 0;
   const Bit Actor::Flags::SetupPostNeighbor = 1;
   const Bit Actor::Flags::SetupPostForce = 2;
   const Bit Actor::Flags::PreIntegrate1 = 3;
   const Bit Actor::Flags::PostIntegrate1 = 4;
   const Bit Actor::Flags::PreTransform = 5;
   const Bit Actor::Flags::PreExchange = 6;
   const Bit Actor::Flags::PostExchange = 7;
   const Bit Actor::Flags::PostNeighbor = 8;
   const Bit Actor::Flags::PreUpdate = 9;
   const Bit Actor::Flags::PostUpdate = 10;
   const Bit Actor::Flags::PreForce = 12;
   const Bit Actor::Flags::PostForce = 13;
   const Bit Actor::Flags::EndOfStep = 14;
   const Bit Actor::Flags::PackExchange = 15;
   const Bit Actor::Flags::PackUpdate = 16;
   const Bit Actor::Flags::PackReverseUpdate = 17;

   /*
   * Constructor.
   */
   Actor::Actor(Simulation& simulation)
    : flags_(0),
      interval_(1),
      simulationPtr_(&simulation)
   {}

   /*
   * Destructor.
   */
   Actor::~Actor()
   {}

   /*
   * Return true if a flag is set, false otherwise.
   */
   bool Actor::isSet(Bit flag)
   {  flag.isSet(flags_); }

   /*
   * Return true if a flag is set, false otherwise.
   */
   void Actor::set(Bit flag)
   {  flag.set(flags_); }

}
#endif
