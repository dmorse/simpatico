#ifndef DDMD_PAIR_FACTORY_H
#define DDMD_PAIR_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>                  // base class template
#include <ddMd/potentials/pair/PairPotential.h>  // template argument

#include <string>
#include <vector>

namespace DdMd
{

   class Simulation;

   /**
   * Factory for PairPotential objects.
   *
   * \ingroup Pair_Module
   */
   class PairFactory : public Factory<PairPotential>
   {

   public:
   
      /**
      * Default constructor.
      */
      PairFactory(Simulation& simulation);

      /**
      * Return a pointer to a new McPairInteration, if possible.
      */
      PairPotential* factory(const std::string& subclass) const;

   private:

      // Pointer to the parent Simulation.
      Simulation* simulationPtr_;

   };
  
}
#endif
