/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesManager.h"
#include "SpeciesFactory.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   SpeciesManager::SpeciesManager()
    : Manager<Species>()
   {  setClassName("SpeciesManager"); }
   
   /*
   * Return a pointer to a new SpeciesFactory object.
   */
   Factory<Simp::Species>* SpeciesManager::newDefaultFactory() const
   { return new SpeciesFactory(); }

}
