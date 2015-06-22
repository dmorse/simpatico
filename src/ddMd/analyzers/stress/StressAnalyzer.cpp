/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StressAnalyzer.h"
#include <ddMd/simulation/Simulation.h>
#include <util/space/Tensor.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   StressAnalyzer::StressAnalyzer(Simulation& simulation) 
    : SymmTensorAverageAnalyzer(simulation)
   {  setClassName("StressAnalyzer"); }

   /*
   * Destructor.
   */
   StressAnalyzer::~StressAnalyzer() 
   {}  

   /*
   * Compute current value.
   */
   void StressAnalyzer::compute() 
   {  
      simulation().computeVirialStress(); 
      simulation().computeKineticStress(); 
   }

   /*
   * Get value current value (call only on master)
   */
   Tensor StressAnalyzer::value() 
   {
      if (!simulation().domain().isMaster()) {
         UTIL_THROW("Error: Not master processor");
      }
      Tensor stress;
      stress = simulation().virialStress();
      stress += simulation().kineticStress();
      return stress;
   }

}
