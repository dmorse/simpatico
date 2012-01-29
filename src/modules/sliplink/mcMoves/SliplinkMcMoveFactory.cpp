#ifndef SLIPLINK_MC_MOVE_FACTORY_CPP
#define SLIPLINK_MC_MOVE_FACTORY_CPP

#include "SliplinkMcMoveFactory.h" 
#include <mcMd/mcSimulation/McSystem.h> 

// Include headers for user defined McMoves
#include "SliplinkerAll.h"
#include "SliplinkerEnd.h"
#include "SliplinkMove.h"
#include "GcSliplinkMove.h"

namespace McMd
{


   /*
   * Constructor.
   */
   SliplinkMcMoveFactory::SliplinkMcMoveFactory(McSimulation& simulation, 
                                                McSystem& system)
    : McMoveFactory(simulation, system)
   {}

   /* 
   * Return a pointer to a new instance of className.
   */
   McMove* SliplinkMcMoveFactory::factory(const std::string &className) const
   {
      McMove* ptr = 0;

      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "SliplinkerAll") {
         ptr = new SliplinkerAll(system());
      } else     
      if (className == "SliplinkerEnd") {
         ptr = new SliplinkerEnd(system());
      } else
      if (className == "SliplinkMove") {
        ptr = new SliplinkMove(system());
      } else          
      if (className == "GcSliplinkMove") {
        ptr = new GcSliplinkMove(system());
      }          
      
      #if 0     
      // If not a user-defined class, try the standard factory 
      if (!ptr) {
         ptr = McMoveFactory::factory(className);
      }
      #endif

      return ptr;
   }

}

#endif
