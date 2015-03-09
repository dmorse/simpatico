/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExternalEnergyAnalyzer.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/potentials/external/ExternalPotential.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   ExternalEnergyAnalyzer::ExternalEnergyAnalyzer(Simulation& simulation) 
    : AverageAnalyzer(simulation)
   {  setClassName("ExternalEnergyAnalyzer"); }

   /*
   * Destructor.
   */
   ExternalEnergyAnalyzer::~ExternalEnergyAnalyzer() 
   {}  

   /*
   * Compute current value.
   */
   void ExternalEnergyAnalyzer::compute() 
   {
      MPI::Intracomm& communicator = simulation().domain().communicator();  
      simulation().externalPotential().computeEnergy(communicator); 
   }

   /*
   * Get value current value (call only on master)
   */
   double ExternalEnergyAnalyzer::value() 
   {
      if (!simulation().domain().isMaster()) {
         UTIL_THROW("Error: Not master processor");
      }
      return simulation().externalPotential().energy();
   }

}
