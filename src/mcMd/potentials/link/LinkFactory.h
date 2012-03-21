#ifdef  MCMD_LINK
#ifndef MCMD_LINK_FACTORY_H
#define MCMD_LINK_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>
#include <iostream>

namespace McMd
{

   using namespace Util;

   class System;
   class BondPotential;

   /**
   * Factory for subclasses of MdBondPotential or McBondPotential.
   * 
   * \ingroup Link_Module
   */
   class LinkFactory : public Factory<BondPotential>
   {
   
   public:
   
      /**
      * Default constructor.
      */
      LinkFactory(System& system);

      /**
      * Return a pointer to a new McBondInteration, if possible.
      */
      BondPotential* factory(const std::string& subclass) const;

   private:

      // Pointer to parent System. 
      System* systemPtr_;

   };
  
}
#endif
#endif
