#ifdef  MCMD_PERTURB
#ifndef HOOMD_MC_PERTURBATION_FACTORY_H
#define HOOMD_MC_PERTURBATION_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
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
   * \ingroup Perturb_Module
   */
   class HoomdMcPerturbationFactory : public Factory<Perturbation> 
   {

   public:

      /**
      * Constructor
      *
      * \param system Parent McSystem
      */
      HoomdMcPerturbationFactory(McSystem& system);

      /**
      * Method to create the hoomd-specific pair perturbations
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
#endif  // #ifdef  MCMD_PERTURB
