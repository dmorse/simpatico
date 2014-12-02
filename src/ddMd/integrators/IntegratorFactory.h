#ifndef DDMD_INTEGRATOR_FACTORY_H
#define DDMD_INTEGRATOR_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <ddMd/integrators/Integrator.h>

#include <string>

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Factory for subclasses of Integrator (i.e., MD integrators).
   *
   * \ingroup DdMd_Integrator_Module
   */
   class IntegratorFactory : public Factory<Integrator> 
   {

   public:

      /// Constructor
      IntegratorFactory(Simulation& simulation);

      /**
      * Method to create any species supplied with Simpatico.
      *
      * \param speciesName name of the Integrator subclass
      * \return Integrator* pointer to new instance of speciesName
      */
      Integrator* factory(const std::string &speciesName) const;

   private:

      /// Pointer to a parent Simulation.
      Simulation* simulationPtr_;

   };

}
#endif
