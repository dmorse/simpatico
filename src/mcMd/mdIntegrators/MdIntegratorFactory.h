#ifndef MCMD_MD_INTEGRATOR_FACTORY_H
#define MCMD_MD_INTEGRATOR_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <mcMd/mdIntegrators/MdIntegrator.h>

#include <string>

namespace McMd
{

   using namespace Util;

   /**
   * Default Factory for subclasses of MdIntegrator.
   *
   * \sa \ref mcMd_mdIntegrators_page
   */
   class MdIntegratorFactory : public Factory<MdIntegrator> 
   {

   public:

      /// Constructor
      MdIntegratorFactory(MdSystem& system);

      /**
      * Method to create any species supplied with Simpatico.
      *
      * \param speciesName name of the MdIntegrator subclass
      * \return MdIntegrator* pointer to new instance of speciesName
      */
      MdIntegrator* factory(const std::string &speciesName) const;

   private:

      /// Pointer to a parent MdSystem.
      MdSystem* systemPtr_;

   };

}
#endif
