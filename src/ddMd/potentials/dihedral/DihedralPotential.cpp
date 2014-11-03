/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
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
