#ifndef MCMD_EXTERNAL_FACADE_H
#define MCMD_EXTERNAL_FACADE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/special/MdPotentialFacade.h>   // base class
#include <mcMd/potentials/external/ExternalPotential.h>  // templ arg
#include <mcMd/potentials/external/ExternalFactory.h>    // templ arg

namespace McMd
{

   class System;

   using namespace Util;

   /**
   * Potential for testing purposes.
   *
   * \ingroup McMd_Potential_Module
   */
   class MdExternalFacade : 
         public MdPotentialFacade<ExternalPotential, ExternalFactory>
   {

   public:

      /**
      * Constructor.
      */
      MdExternalFacade(System& system);

   };

} 
#endif
