#ifndef DDMD_INTEGRATOR_FACTORY_CPP
#define DDMD_INTEGRATOR_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "IntegratorFactory.h"  

// Subclasses of Integrator 
#include "NveIntegrator.h"
#include "NvtIntegrator.h"
//#include "NvtDpdVvIntegrator.h"
//#include "NphIntegrator.h"

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor
   */
   IntegratorFactory::IntegratorFactory(System& system)
    : systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Integrator subclass className.
   */
   Integrator* IntegratorFactory::factory(const std::string &className) const
   {
      Integrator *ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "NveIntegrator") {
         ptr = new NveIntegrator(*systemPtr_);
      } else
      if (className == "NvtIntegrator") {
         ptr = new NvtIntegrator(*systemPtr_);
      } // else
      //if (className == "NvtDpdVvIntegrator") {
      //   ptr = new NvtDpdVvIntegrator(*systemPtr_);
      //}
      //if (className == "NphIntegrator") {
      //   ptr = new NphIntegrator(*systemPtr_);
      //}
      return ptr;
   }

}
#endif
