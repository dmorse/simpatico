#ifndef MCMD_DIHEDRAL_FACTORY_H
#define MCMD_DIHEDRAL_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>
#include <iostream>

namespace McMd
{

   using namespace Util;

   class System;
   class DihedralPotential;

   /**
   * Factory for subclasses of DihedralPotential.
   * 
   * \ingroup McMd_Dihedral_Module
   */
   class DihedralFactory : public Factory<DihedralPotential>
   {
   
   public:
   
      /**
      * Default constructor.
      */
      DihedralFactory(System& system);

      /**
      * Return a pointer to a new DihedralPotential, if possible.
      */
      DihedralPotential* factory(const std::string& subclass) const;

   private:

      // Pointer to parent System. 
      System* systemPtr_;

   };
  
}

#endif // ifndef DIHEDRAL_FACTORY_H
