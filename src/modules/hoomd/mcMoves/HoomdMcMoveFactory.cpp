#ifndef HOOMD_MC_MOVE_FACTORY_CPP
#define HOOMD_MC_MOVE_FACTORY_CPP

#include "HoomdMcMoveFactory.h" 
#include <mcMd/mcSimulation/McSystem.h> 

#include <modules/hoomd/mcMoves/HoomdMove.h>
#include <modules/hoomd/mcMoves/HoomdNPTMTKMove.h>

namespace McMd
{


   /*
   * Constructor.
   */
   HoomdMcMoveFactory::HoomdMcMoveFactory(McSimulation& simulation, 
                                                McSystem& system)
    : McMoveFactory(simulation, system)
   {}

   /* 
   * Return a pointer to a new instance of className.
   */
   McMove* HoomdMcMoveFactory::factory(const std::string &className) const
   {
      McMove* spp = 0;

      if (className == "HoomdMove") {
        spp = new HoomdMove(system());
      }           
//      if (className == "HoomdNPTMTKMove") {  //Temporarily disabled
//        spp = new HoomdNPTMTKMove(system());           
//      }
      return spp;
   }

}

#endif
