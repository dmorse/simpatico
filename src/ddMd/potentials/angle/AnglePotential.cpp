#ifndef DDMD_ANGLE_POTENTIAL_CPP
#define DDMD_ANGLE_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AnglePotential.h"
#include <ddMd/simulation/Simulation.h>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   AnglePotential::AnglePotential(Simulation& simulation)
    : boundaryPtr_(&simulation.boundary()),
      storagePtr_(&simulation.angleStorage())
   { setClassName("AnglePotential"); }

   /*
   * Default constructor (for unit testing).
   */
   AnglePotential::AnglePotential()
    : boundaryPtr_(0),
      storagePtr_(0)
   { setClassName("AnglePotential"); }

   /*
   * Associate with related objects. (for unit testing).
   */
   void AnglePotential::associate(Boundary& boundary, GroupStorage<3>& storage)
   {
      boundaryPtr_ = &boundary;
      storagePtr_ = &storage;
   } 

   /*
   * Destructor.
   */
   AnglePotential::~AnglePotential()
   {}

}
#endif
