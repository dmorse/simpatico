#ifndef DIAGNOSTIC_FACTORY
#define DIAGNOSTIC_FACTORY

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DiagnosticFactory.h" // Class header

//#include <ddMd/simulation/Simulation.h>

// Diagnostics 
#include "WriteConfig.h"
#include "OutputEnergy.h"
#include "OutputPressure.h"
#include "OutputBoxdim.h"
#include "OutputTemperature.h"
#include "OutputPairEnergies.h"
#include "StructureFactor.h"
#include "StructureFactorGrid.h"

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DiagnosticFactory::DiagnosticFactory(Simulation& simulation)
    : simulationPtr_(&simulation)
   {}

   /* 
   * Return a pointer to an instance of Diagnostic subclass className.
   */
   Diagnostic* DiagnosticFactory::factory(const std::string &className) const
   {
      Diagnostic* ptr = 0;

      // Try subfactories first (if any)
      //ptr = trySubfactories(className);
      //if (ptr) return ptr;

      // Simulation Diagnostics
      if (className == "WriteConfig") {
         ptr = new WriteConfig(simulation());
      } else 
      if (className == "OutputEnergy") {
         ptr = new OutputEnergy(simulation());
      } else 
      if (className == "OutputPressure") {
         ptr = new OutputPressure(simulation());
      } else
      if (className == "OutputBoxdim") {
         ptr = new OutputBoxdim(simulation());
      } else
      if (className == "OutputTemperature") {
         ptr = new OutputTemperature(simulation());
      } else
      if (className == "OutputPairEnergies") {
         ptr = new OutputPairEnergies(simulation());
      } else
      if (className == "StructureFactor") {
         ptr = new StructureFactor(simulation());
      } else
      if (className == "StructureFactorGrid") {
         ptr = new StructureFactorGrid(simulation());
      }
      return ptr;
   }

}
#endif
