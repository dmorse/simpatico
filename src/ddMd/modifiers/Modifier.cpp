#ifndef DDMD_ACTOR_CPP
#define DDMD_ACTOR_CPP

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

   const Bit Modifier::Flags::SetupPostExchange = 0;
   const Bit Modifier::Flags::SetupPostNeighbor = 1;
   const Bit Modifier::Flags::SetupPostForce = 2;
   const Bit Modifier::Flags::PreIntegrate1 = 3;
   const Bit Modifier::Flags::PostIntegrate1 = 4;
   const Bit Modifier::Flags::PreTransform = 5;
   const Bit Modifier::Flags::PreExchange = 6;
   const Bit Modifier::Flags::PostExchange = 7;
   const Bit Modifier::Flags::PostNeighbor = 8;
   const Bit Modifier::Flags::PreUpdate = 9;
   const Bit Modifier::Flags::PostUpdate = 10;
   const Bit Modifier::Flags::PreForce = 12;
   const Bit Modifier::Flags::PostForce = 13;
   const Bit Modifier::Flags::EndOfStep = 14;
   const Bit Modifier::Flags::Exchange = 15;
   const Bit Modifier::Flags::Update = 16;
   const Bit Modifier::Flags::ReverseUpdate = 17;

   /*
   * Default constructor.
   */
   Modifier::Modifier()
    : flags_(0),
      interval_(1),
      simulationPtr_(0)
   {}

   /*
   * Constructor.
   */
   Modifier::Modifier(Simulation& simulation)
    : flags_(0),
      interval_(1),
      simulationPtr_(&simulation)
   {}

   /*
   * Destructor.
   */
   Modifier::~Modifier()
   {}

   /*
   * Return true if a flag is set, false otherwise.
   */
   bool Modifier::isSet(Bit flag)
   {  return flag.isSet(flags_); }

   /*
   * Return true if a flag is set, false otherwise.
   */
   void Modifier::set(Bit flag)
   {  flag.set(flags_); }

}
#endif
