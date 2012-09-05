#ifndef MCMD_SPECIES_MANAGER_H
#define MCMD_SPECIES_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>
#include <util/param/Manager.h>
#include "Species.h"

namespace McMd
{

   using namespace Util;

   /**
   * A Manager for a set of Species objects.
   *
   * \ingroup McMd_Manager_Module
   * \ingroup McMd_Species_Module
   */
   class SpeciesManager : public Manager<Species>
   {

   public:

      /**
      * Constructor.
      */
      SpeciesManager();

   protected:

      virtual Factory<Species>* newDefaultFactory() const;

   };

}
#endif
