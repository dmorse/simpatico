#ifndef SPECIES_MANAGER_CPP
#define SPECIES_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesManager.h"
#include "SpeciesFactory.h"

namespace McMd
{

   using namespace Util;

   /*
   * Read parameter file. 
   */
   void SpeciesManager::readParam(std::istream &in)
   {
      readBegin(in, "SpeciesManager");
      Manager<Species>::readParam(in);
   } 

   /*
   * Return a pointer to a new SpeciesFactory object.
   */
   Factory<Species>* SpeciesManager::newDefaultFactory() const
   { return new SpeciesFactory(); }

}

#endif
