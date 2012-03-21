#ifndef DIAGNOSTIC_FACTORY
#define DIAGNOSTIC_FACTORY

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DiagnosticFactory.h" // Class header

//#include <ddMd/system/System.h>

// Diagnostics for MdSystem only
#include "WriteConfig.h"
#include "OutputEnergy.h"

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DiagnosticFactory::DiagnosticFactory(System& system)
    : systemPtr_(&system)
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

      // System Diagnostics
      if (className == "WriteConfig") {
         ptr = new WriteConfig(system());
      } else 
      if (className == "OutputEnergy") {
         ptr = new OutputEnergy(system());
      }
      return ptr;
   }

}
#endif
