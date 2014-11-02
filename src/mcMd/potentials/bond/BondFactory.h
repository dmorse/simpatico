#ifndef MCMD_BOND_FACTORY_H
#define MCMD_BOND_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
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
   * \ingroup McMd_Bond_Module
   */
   class BondFactory : public Factory<BondPotential>
   {

   public:
   
      /**
      * Default constructor.
      */
      BondFactory(System& system);

      /**
      * Return a pointer to a new BondPotential, if possible.
      */
      BondPotential* factory(const std::string& subclass) const;

   private:

      // Pointer to the parent System.
      System* systemPtr_;

   };
  
}
#endif
