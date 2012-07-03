#ifndef MCMD_MD_INTEGRATOR_FACTORY_CPP
#define MCMD_MD_INTEGRATOR_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdIntegratorFactory.h"  

// Subclasses of MdIntegrator 
#include "NveVvIntegrator.h"
#include "NvtNhIntegrator.h"
#include "NvtDpdVvIntegrator.h"
#include "NphIntegrator.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor
   */
   MdIntegratorFactory::MdIntegratorFactory(MdSystem& system)
    : systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of MdIntegrator subclass className.
   */
   MdIntegrator* MdIntegratorFactory::factory(const std::string &className) const
   {
      MdIntegrator *ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "NveVvIntegrator") {
         ptr = new NveVvIntegrator(*systemPtr_);

      } else
      if (className == "NvtNhIntegrator") {
         ptr = new NvtNhIntegrator(*systemPtr_);
      } else
      if (className == "NvtDpdVvIntegrator") {
         ptr = new NvtDpdVvIntegrator(*systemPtr_);
      }
      if (className == "NphIntegrator") {
         ptr = new NphIntegrator(*systemPtr_);
      }
      return ptr;
   }

}
#endif
