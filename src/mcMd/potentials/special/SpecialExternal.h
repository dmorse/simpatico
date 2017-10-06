#ifndef MCMD_SPECIAL_EXTERNAL_H
#define MCMD_SPECIAL_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/special/SpecialPotentialFacade.h>   // base class
#include <mcMd/potentials/external/ExternalPotential.h>       // templ arg
#include <mcMd/potentials/external/ExternalFactory.h>         // templ arg

namespace McMd
{

   class System;

   using namespace Util;

   /**
   * Potential for testing purposes.
   *
   * \ingroup McMd_Potential_Module
   */
   class SpecialExternal : 
         public SpecialPotentialFacade<ExternalPotential, ExternalFactory>
   {

   public:

      /**
      * Constructor.
      */
      SpecialExternal(System& system);

   };

} 
#endif
