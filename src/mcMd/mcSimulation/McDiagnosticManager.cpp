#ifndef MCMD_MC_DIAGNOSTIC_MANAGER_CPP
#define MCMD_MC_DIAGNOSTIC_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McDiagnosticManager.h"                
#include "McSimulation.h"                        
#include <mcMd/diagnostics/mcSystem/McDiagnosticFactory.h> 

namespace McMd
{

   using namespace Util;

   // Constructor.
   McDiagnosticManager::McDiagnosticManager(McSimulation& simulation)
    : simulationPtr_(&simulation),
      systemPtr_(&simulation.system())
   {}

   // Constructor.
   McDiagnosticManager::McDiagnosticManager(McSimulation& simulation, 
		                            McSystem &system)
    : simulationPtr_(&simulation),
      systemPtr_(&system)
   {}

   /// Return pointer to a new DiagnosticFactory.
   Factory<Diagnostic>* McDiagnosticManager::newDefaultFactory() const
   {  return new McDiagnosticFactory(*simulationPtr_, *systemPtr_); }

}
#endif
