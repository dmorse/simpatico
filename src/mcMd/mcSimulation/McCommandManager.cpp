/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McCommandManager.h"                
#include "McSimulation.h"  
#include <mcMd/commands/McCommandFactory.h> 

namespace McMd
{

   using namespace Util;

   // Constructor.
   McCommandManager::McCommandManager(McSimulation& simulation)
    : simulationPtr_(&simulation),
      systemPtr_(&simulation.system())
   {  setClassName("McCommandManager"); }

   // Constructor.
   McCommandManager::McCommandManager(McSimulation& simulation, 
		                            McSystem& system)
    : simulationPtr_(&simulation),
      systemPtr_(&system)
   {}

   // Destructor.
   McCommandManager::~McCommandManager()
   {}

   /// Return pointer to a new CommandFactory.
   Factory<Command>* McCommandManager::newDefaultFactory() const
   {  return new McCommandFactory(*simulationPtr_, *systemPtr_); }

}
