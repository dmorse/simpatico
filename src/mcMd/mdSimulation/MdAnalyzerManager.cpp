/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdAnalyzerManager.h"                
#include "MdSimulation.h"  
#include <mcMd/analyzers/mdSystem/MdAnalyzerFactory.h> 

namespace McMd
{

   using namespace Util;

   // Constructor.
   MdAnalyzerManager::MdAnalyzerManager(MdSimulation& simulation)
    : simulationPtr_(&simulation),
      systemPtr_(&simulation.system())
   {}

   // Constructor.
   MdAnalyzerManager::MdAnalyzerManager(MdSimulation& simulation, 
		                            MdSystem& system)
    : simulationPtr_(&simulation),
      systemPtr_(&system)
   {}

   // Destructor.
   MdAnalyzerManager::~MdAnalyzerManager()
   {}

   /// Return pointer to a new AnalyzerFactory.
   Factory<Analyzer>* MdAnalyzerManager::newDefaultFactory() const
   {  return new MdAnalyzerFactory(*simulationPtr_, *systemPtr_); }

}
