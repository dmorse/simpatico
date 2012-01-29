#ifndef MC_MODULE_CPP
#define MC_MODULE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McModule.h"
#include <mcMd/mcSimulation/McSimulation.h>

namespace McMd
{

   /*
   * Constructor.
   */
   McModule::McModule(McSimulation& simulation) 
    : simulationPtr_(&simulation),
      systemPtr_(&simulation.system())
   {}

   /*
   * Destructor.
   */
   McModule::~McModule()
   {}

   /*
   * Return McSimulation by reference.
   */
   McSimulation& McModule::simulation() const
   {  return *simulationPtr_; }

   /*
   * Return McSystem by reference.
   */
   McSystem& McModule::system() const
   {  return *systemPtr_; }

}
#endif
