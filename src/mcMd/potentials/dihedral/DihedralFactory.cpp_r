#ifndef MCMD_DIHEDRAL_FACTORY_CPP
#define MCMD_DIHEDRAL_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/dihedral/DihedralFactory.h>
#include <mcMd/simulation/System.h>
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#include <mcMd/potentials/dihedral/DihedralPotentialImpl.h>

// Dihedral interaction classes
#include <inter/dihedral/CosineDihedral.h>

namespace McMd
{

   using namespace Inter;

   /**
   * Default constructor.
   */
   DihedralFactory::DihedralFactory(System& system)
    : systemPtr_(&system)
   {}

   /*
   * Return a pointer to a new DihedralPotential, if possible.
   */
   DihedralPotential* DihedralFactory::factory(const std::string& name) const
   {
      DihedralPotential* ptr = 0;

      ptr = trySubfactories(name);
      if (ptr) return ptr;

      if (name == "CosineDihedral") {
         ptr = new DihedralPotentialImpl<CosineDihedral>(*systemPtr_);
      }
      return ptr;
   }

}
#endif // ifndef DIHEDRAL_FACTORY_CPP
