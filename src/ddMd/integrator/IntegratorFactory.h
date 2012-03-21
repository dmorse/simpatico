#ifndef DDMD_INTEGRATOR_FACTORY_H
#define DDMD_INTEGRATOR_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <ddMd/integrator/Integrator.h>

#include <string>

namespace DdMd
{

   class System;

   using namespace Util;

   /**
   * Default Factory for subclasses of Integrator.
   */
   class IntegratorFactory : public Factory<Integrator> 
   {

   public:

      /// Constructor
      IntegratorFactory(System& system);

      /**
      * Method to create any species supplied with Simpatico.
      *
      * \param speciesName name of the Integrator subclass
      * \return Integrator* pointer to new instance of speciesName
      */
      Integrator* factory(const std::string &speciesName) const;

   private:

      /// Pointer to a parent System.
      System* systemPtr_;

   };

}
#endif
