#ifndef SPECIES_MANAGER_H
#define SPECIES_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
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
   * \ingroup Manager_Module
   * \ingroup Species_Module
   */
   class SpeciesManager : public Manager<Species>
   {

   public:

      /// Constructor.
      SpeciesManager()
       : Manager<Species>()
      {}

      /**
      * Read parameter file. 
      *
      * \param in input parameter file stream.
      */
      virtual void readParam(std::istream &in);

   protected:

      virtual Factory<Species>* newDefaultFactory() const;

   };

}
#endif
