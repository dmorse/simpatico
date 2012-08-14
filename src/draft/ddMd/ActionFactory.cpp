#ifndef ACTION_FACTORY
#define ACTION_FACTORY

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ActionFactory.h" // Class header

// Actions 
// #include "WriteConfig.h"
// #include "OutputEnergy.h"

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   ActionFactory::ActionFactory(Simulation& simulation)
    : simulationPtr_(&simulation)
   {}

   /* 
   * Return a pointer to an instance of Action subclass className.
   */
   Action* ActionFactory::factory(const std::string &className) const
   {
      Action* ptr = 0;

      // Try subfactories first (if any)
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      #if 0
      // Simulation Actions
      if (className == "WriteConfig") {
         ptr = new WriteConfig(simulation());
      } // else 
      #endif

      return ptr;
   }

}
#endif
