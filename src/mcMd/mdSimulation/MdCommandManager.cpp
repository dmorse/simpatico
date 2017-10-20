/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdCommandManager.h"                
#include "MdSimulation.h"  
#include <mcMd/commands/MdCommandFactory.h> 

namespace McMd
{

   using namespace Util;

   // Constructor.
   MdCommandManager::MdCommandManager(MdSimulation& simulation)
    : CommandManager(),
      simulationPtr_(&simulation),
      systemPtr_(&simulation.system())
   {  setClassName("MdCommandManager"); }

   // Constructor.
   MdCommandManager::MdCommandManager(MdSimulation& simulation, 
		                            MdSystem& system)
    : CommandManager(),
      simulationPtr_(&simulation),
      systemPtr_(&system)
   {}

   // Destructor.
   MdCommandManager::~MdCommandManager()
   {}

   /// Return pointer to a new CommandFactory.
   Factory<Command>* MdCommandManager::newDefaultFactory() const
   {  return new MdCommandFactory(*simulationPtr_, *systemPtr_); }

}
