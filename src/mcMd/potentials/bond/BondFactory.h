#ifndef BOND_FACTORY_H
#define BOND_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>                  // base class template
#include <mcMd/potentials/bond/BondPotential.h>  // template argument

#include <string>
#include <vector>

namespace McMd
{

   class System;

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
      BondFactory(System& system);

      /**
      * Return a pointer to a new McBondInteration, if possible.
      */
      BondPotential* factory(const std::string& subclass) const;

   private:

      // Pointer to the parent System.
      System* systemPtr_;

   };
  
}
#endif
