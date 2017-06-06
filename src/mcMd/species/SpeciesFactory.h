#ifndef MCMD_SPECIES_FACTORY_H
#define MCMD_SPECIES_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <simp/species/Species.h>

#include <string>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * Default Factory for subclasses of Species.
   */
   class SpeciesFactory : public Factory<Simp::Species> 
   {

   public:

      /**
      * Method to create any species supplied with Simpatico.
      *
      * \param speciesName name of the Species subclass
      * \return Species* pointer to new instance of speciesName
      */
      Simp::Species* factory(const std::string &speciesName) const;

   };

}
#endif
