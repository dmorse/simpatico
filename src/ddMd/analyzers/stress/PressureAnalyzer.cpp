/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PressureAnalyzer.h"
#include <ddMd/simulation/Simulation.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   PressureAnalyzer::PressureAnalyzer(Simulation& simulation) 
    : AverageAnalyzer(simulation)
   {  setClassName("PressureAnalyzer"); }

   /*
   * Destructor.
   */
   PressureAnalyzer::~PressureAnalyzer() 
   {}  

   /*
   * Compute current value.
   */
   void PressureAnalyzer::compute() 
   {  
      simulation().computeVirialStress(); 
      simulation().computeKineticStress(); 
   }

   /*
   * Get value current value (call only on master)
   */
   double PressureAnalyzer::value() 
   {
      if (!simulation().domain().isMaster()) {
         UTIL_THROW("Error: Not master processor");
      }
      double pressure;
      pressure = simulation().virialPressure();
      pressure += simulation().kineticPressure();
      return pressure;
   }

}
