#ifndef DDMD_DIHEDRAL_POTENTIAL_CPP
#define DDMD_DIHEDRAL_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DihedralPotential.h"
#include <ddMd/simulation/Simulation.h>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   DihedralPotential::DihedralPotential(Simulation& simulation)
    : boundaryPtr_(&simulation.boundary()),
      storagePtr_(&simulation.dihedralStorage())
   {  setClassName("DihedralPotential");  }

   /*
   * Default constructor (for unit testing).
   */
   DihedralPotential::DihedralPotential()
    : boundaryPtr_(0),
      storagePtr_(0)
   {  setClassName("DihedralPotential");  }

   /*
   * Associate with related objects (for unit testing).
   */
   void DihedralPotential::associate(Boundary& boundary, GroupStorage<4>& storage)
   {
      boundaryPtr_ = &boundary;
      storagePtr_ = &storage;
   } 

   /*
   * Destructor.
   */
   DihedralPotential::~DihedralPotential()
   {}

}
#endif
