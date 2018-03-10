#ifndef MCMD_MC_PERTURBATION_FACTORY_H
#define MCMD_MC_PERTURBATION_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <mcMd/perturb/Perturbation.h>

#include <string>

namespace McMd
{

   using namespace Util;
   class McSystem;

   /**
   * Default Factory for subclasses of Perturbation.
   *
   * \ingroup McMd_Perturb_Module
   */
   class McPerturbationFactory : public Factory<Perturbation> 
   {

   public:

      /**
      * Constructor
      *
      * \param system Parent McSystem
      */
      McPerturbationFactory(McSystem& system);

      /**
      * Method to create any species supplied with Simpatico.
      *
      * \param  className     name of the Perturbation subclass
      * \return Perturbation* pointer to new instance of className
      */
      Perturbation* factory(const std::string &className) const;

   private:

      McSystem* systemPtr_;

   };

}

#endif
