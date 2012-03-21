#ifndef MCMD_SPECIES_FACTORY_H
#define MCMD_SPECIES_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include "Species.h"

#include <string>

namespace McMd
{

   using namespace Util;

   /**
   * Default Factory for subclasses of Species.
   */
   class SpeciesFactory : public Factory<Species> 
   {

   public:

      /**
      * Method to create any species supplied with Simpatico.
      *
      * \param speciesName name of the Species subclass
      * \return Species* pointer to new instance of speciesName
      */
      Species* factory(const std::string &speciesName) const;

   };

}
#endif
