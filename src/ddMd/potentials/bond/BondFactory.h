#ifndef DDMD_BOND_FACTORY_H
#define DDMD_BOND_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>                  // base class template
#include <ddMd/potentials/bond/BondPotential.h>  // template argument

#include <string>

namespace DdMd
{

   class Simulation;

   /**
   * Factory for BondPotential objects.
   *
   * \ingroup Bond_Module
   */
   class BondFactory : public Factory<BondPotential>
   {

   public:
   
      /**
      * Default constructor.
      */
      BondFactory(Simulation& simulation);

      /**
      * Return a pointer to a new McBondInteration, if possible.
      */
      BondPotential* factory(const std::string& subclass) const;

   private:

      // Pointer to the parent Simulation.
      Simulation* simulationPtr_;

   };
  
}
#endif
