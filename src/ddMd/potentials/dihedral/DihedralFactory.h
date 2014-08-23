#ifndef DDMD_DIHEDRAL_FACTORY_H
#define DDMD_DIHEDRAL_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>                         // base class template
#include <ddMd/potentials/dihedral/DihedralPotential.h> // template argument

#include <string>

namespace DdMd
{

   class Simulation;

   /**
   * Factory for DihedralPotential objects.
   *
   * \ingroup DdMd_Dihedral_Module
   */
   class DihedralFactory : public Factory<DihedralPotential>
   {

   public:
   
      /**
      * Default constructor.
      */
      DihedralFactory(Simulation& simulation);

      /**
      * Return a pointer to a new McDihedralInteration, if possible.
      */
      DihedralPotential* factory(const std::string& subclass) const;

   private:

      // Pointer to the parent Simulation.
      Simulation* simulationPtr_;

   };
  
}
#endif
