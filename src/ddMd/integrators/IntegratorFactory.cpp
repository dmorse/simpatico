#ifndef DDMD_INTEGRATOR_FACTORY_CPP
#define DDMD_INTEGRATOR_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "IntegratorFactory.h"  

// Subclasses of Integrator 
#include "NveIntegrator.h"
#include "NvtIntegrator.h"
#include "NptIntegrator.h"
#include "NphIntegrator.h"
//#include "NvtDpdVvIntegrator.h"
//#include "NphIntegrator.h"

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor
   */
   IntegratorFactory::IntegratorFactory(Simulation& simulation)
    : simulationPtr_(&simulation)
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
         ptr = new NveIntegrator(*simulationPtr_);
      } else
      if (className == "NvtIntegrator") {
         ptr = new NvtIntegrator(*simulationPtr_);
      }
      if (className == "NptIntegrator") {
         ptr = new NptIntegrator(*simulationPtr_);
      }
      if (className == "NphIntegrator") {
         ptr = new NphIntegrator(*simulationPtr_);
      }
      // else
      //if (className == "NvtDpdVvIntegrator") {
      //   ptr = new NvtDpdVvIntegrator(*simulationPtr_);
      //}
      //if (className == "NphIntegrator") {
      //   ptr = new NphIntegrator(*simulationPtr_);
      //}
      return ptr;
   }

}
#endif
