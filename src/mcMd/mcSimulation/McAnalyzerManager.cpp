/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McAnalyzerManager.h"                
#include "McSimulation.h"                        
#include <mcMd/analyzers/mcSystem/McAnalyzerFactory.h> 

namespace McMd
{

   using namespace Util;

   // Constructor.
   McAnalyzerManager::McAnalyzerManager(McSimulation& simulation)
    : simulationPtr_(&simulation),
      systemPtr_(&simulation.system())
   {}

   // Constructor.
   McAnalyzerManager::McAnalyzerManager(McSimulation& simulation, 
		                            McSystem &system)
    : simulationPtr_(&simulation),
      systemPtr_(&system)
   {}

   /// Return pointer to a new AnalyzerFactory.
   Factory<Analyzer>* McAnalyzerManager::newDefaultFactory() const
   {  return new McAnalyzerFactory(*simulationPtr_, *systemPtr_); }

}
