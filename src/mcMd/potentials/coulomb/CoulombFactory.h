#ifndef MCMD_COULOMB_FACTORY_H
#define MCMD_COULOMB_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>                        // base class template
#include <mcMd/potentials/coulomb/MdCoulombPotential.h>  // template argument
#include <mcMd/potentials/coulomb/MdEwaldPotential.h>

#include <string>
#include <vector>

namespace McMd
{

   class System;

   /**
   * Factory for CoulombPotential objects.
   *
   * \ingroup McMd_Bond_Module
   */
   class CoulombFactory : public Factory<MdCoulombPotential>
   {

   public:
   
      /**
      * Default constructor.
      */
      CoulombFactory(System& system);

      /**
      * Return a pointer to a new CoulombPotential, if possible.
      */
      MdCoulombPotential* factory(const std::string& subclass) const;

   private:

      // Pointer to the parent System.
      System* systemPtr_;

   };
  
}
#endif
