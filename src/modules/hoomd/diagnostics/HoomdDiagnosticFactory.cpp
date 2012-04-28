/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdDiagnosticFactory.h" // Class header

#include "GPUStructureFactorGrid.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   HoomdDiagnosticFactory::HoomdDiagnosticFactory(Simulation& simulation, 
                                            System& system)
    : simulationPtr_(&simulation),
      systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Species subclass className.
   */
   Diagnostic* HoomdDiagnosticFactory::factory(const std::string &className) const
   {
      Diagnostic* ptr = 0;

      if (className == "GPUStructureFactorGrid") {
         ptr = new GPUStructureFactorGrid(*systemPtr_);
      } 

      return ptr;
   }

}

