#ifndef DDMD_ACTOR_FACTORY_CPP
#define DDMD_ACTOR_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ActorFactory.h" // Class header

// Actors 
// #include "WriteConfig.h"
// #include "OutputEnergy.h"

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   ActorFactory::ActorFactory(Simulation& simulation)
    : simulationPtr_(&simulation)
   {}

   /* 
   * Return a pointer to an instance of Actor subclass className.
   */
   Actor* ActorFactory::factory(const std::string &className) const
   {
      Actor* ptr = 0;

      // Try subfactories first (if any)
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      #if 0
      // Simulation Actors
      if (className == "WriteConfig") {
         ptr = new WriteConfig(simulation());
      } // else 
      #endif

      return ptr;
   }

}
#endif
