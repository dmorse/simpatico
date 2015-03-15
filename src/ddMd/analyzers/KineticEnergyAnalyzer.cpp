/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "KineticEnergyAnalyzer.h"
#include <ddMd/simulation/Simulation.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   KineticEnergyAnalyzer::KineticEnergyAnalyzer(Simulation& simulation) 
    : AverageAnalyzer(simulation)
   {  setClassName("KineticEnergyAnalyzer"); }

   /*
   * Destructor.
   */
   KineticEnergyAnalyzer::~KineticEnergyAnalyzer() 
   {}  

   /*
   * Compute current value.
   */
   void KineticEnergyAnalyzer::compute() 
   {  simulation().computeKineticEnergy(); }

   /*
   * Get value current value (call only on master)
   */
   double KineticEnergyAnalyzer::value() 
   {
      if (!simulation().domain().isMaster()) {
         UTIL_THROW("Error: Not master processor");
      }
      return simulation().kineticEnergy();
   }

}
