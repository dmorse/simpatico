/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdAnalyzerFactory.h" // Class header

#include "GPUStructureFactorGrid.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   HoomdAnalyzerFactory::HoomdAnalyzerFactory(Simulation& simulation, 
                                            System& system)
    : simulationPtr_(&simulation),
      systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Species subclass className.
   */
   Analyzer* HoomdAnalyzerFactory::factory(const std::string &className) const
   {
      Analyzer* ptr = 0;

      if (className == "GPUStructureFactorGrid") {
         ptr = new GPUStructureFactorGrid(*systemPtr_);
      } 

      return ptr;
   }

}

