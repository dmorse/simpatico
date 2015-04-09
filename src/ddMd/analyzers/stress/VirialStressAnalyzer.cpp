/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VirialStressAnalyzer.h"
#include <ddMd/simulation/Simulation.h>
#include <util/space/Tensor.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   VirialStressAnalyzer::VirialStressAnalyzer(Simulation& simulation) 
    : SymmTensorAverageAnalyzer(simulation)
   {  setClassName("VirialStressAnalyzer"); }

   /*
   * Destructor.
   */
   VirialStressAnalyzer::~VirialStressAnalyzer() 
   {}  

   /*
   * Compute current value.
   */
   void VirialStressAnalyzer::compute() 
   {  
      simulation().computeVirialStress(); 
   }

   /*
   * Get value current value (call only on master)
   */
   Tensor VirialStressAnalyzer::value() 
   {
      if (!simulation().domain().isMaster()) {
         UTIL_THROW("Error: Not master processor");
      }
      return simulation().virialStress();
   }

}
