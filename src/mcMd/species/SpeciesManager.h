#ifndef MCMD_SPECIES_MANAGER_H
#define MCMD_SPECIES_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>
#include <util/param/Manager.h>
#include <simp/species/Species.h>

namespace McMd
{

   using namespace Util;
   using namespace McMd;

   /**
   * A Manager for a set of Species objects.
   *
   * \ingroup Simp_Manager_Module
   * \ingroup Simp_Species_Module
   */
   class SpeciesManager : public Manager<Simp::Species>
   {

   public:

      /**
      * Constructor.
      */
      SpeciesManager();

   protected:

      virtual Factory<Simp::Species>* newDefaultFactory() const;

   };

}
#endif
