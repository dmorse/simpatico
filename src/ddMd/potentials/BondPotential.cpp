#ifndef BOND_POTENTIAL_CPP
#define BOND_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BondPotential.h"
#include <ddMd/system/System.h>
#include <ddMd/storage/GroupStorage.h>
#include <ddMd/storage/GroupIterator.h>
#include <ddMd/boundary/Boundary.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   BondPotential::BondPotential()
    : boundaryPtr_(0),
      interactionPtr_(0),
      storagePtr_(0)
   {}

   /*
   * Constructor.
   */
   BondPotential::BondPotential(System& system)
    : boundaryPtr_(&system.boundary()),
      interactionPtr_(&system.bondInteraction()),
      storagePtr_(&system.bondStorage())
   {}

   /*
   * Destructor.
   */
   BondPotential::~BondPotential()
   {}

   /*
   * Retain pointers to associated objects.
   */
   void BondPotential::associate(Boundary& boundary,
                                 const BondInteraction& interaction,
                                 GroupStorage<2>& storage)
   {
      boundaryPtr_ = &boundary;
      interactionPtr_ = &interaction;
      storagePtr_ = &storage;
   }

   /*
   * Increment atomic forces, without calculating energy.
   */
   void BondPotential::addForces()
   {  addForces(true, false);  }

   /*
   * Increment atomic forces and compute pair energy for this processor.
   */
   void BondPotential::addForces(double& energy)
   {  energy = addForces(true, true);  }

   /*
   * Calculate and return total bond energy for this processor.
   */
   double BondPotential::energy()
   {  return addForces(false, true);  }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   double BondPotential::addForces(bool needForce, bool needEnergy)
   {
      // Preconditions
      //if (!storagePtr_->isInitialized()) {
      //   UTIL_THROW("AtomStorage must be initialized");
      //}

      Vector f;
      double rsq;
      double energy = 0.0;
      GroupIterator<2> iter;
      Atom* atom0Ptr;
      Atom* atom1Ptr;
      int type;

      storagePtr_->begin(iter);
      for ( ; !iter.atEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         // Set f = r0 - r1, minimum image separation between atoms
         rsq = boundaryPtr_->distanceSq(atom0Ptr->position(), 
                                        atom1Ptr->position(), f);
         if (needEnergy) {
            energy += interactionPtr_->energy(rsq, type);
         }
         if (needForce) {
            // Set force = (r0-r1)*(forceOverR)
            f *= interactionPtr_->forceOverR(rsq, type);
            atom0Ptr->force() += f;
            atom1Ptr->force() -= f;
         }
      }
      return energy;
   }

}
#endif
