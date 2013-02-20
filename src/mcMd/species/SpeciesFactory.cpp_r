#ifndef MCMD_SPECIES_FACTORY_CPP
#define MCMD_SPECIES_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesFactory.h"  

// Subclasses of Species 
#include "Point.h"
#include "Homopolymer.h"
#include "Diblock.h"
#include "HomoRing.h"
#include "HomopolymerSG.h"

namespace McMd
{

   using namespace Util;

   /* 
   * Return a pointer to a instance of Species subclass speciesName.
   */
   Species* SpeciesFactory::factory(const std::string &className) const
   {
      Species *ptr = 0;

      ptr = trySubfactories(className);
      if (ptr) return ptr;     

      // Try explicit class names
      if (className == "Species") {
         ptr = new Species();
      } else
      if (className == "Point") {
         ptr = new Point();
      } else
      if (className == "Homopolymer") {
         ptr = new Homopolymer();
      } else
      if (className == "Diblock") {
         ptr = new Diblock();
      } else
      if (className == "HomoRing") {
         ptr = new HomoRing();
      } else 
      if (className == "HomopolymerSG") {
         ptr = new HomopolymerSG();
      } 
      return ptr;
   }

}
#endif
