#ifndef BOND_POTENTIAL_H
#define BOND_POTENTIAL_H

//#include <util/param/ParamComposite.h>      // base class
#include <ddMd/boundary/Boundary.h>           // typedef
#include <ddMd/potentials/BondInteraction.h>  // typedef

#include <iostream>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   using namespace Util;

   class System;
   template <int N> class GroupStorage;

   /**
   * An BondPotential calculates bond forces and energies for a parent System.
   *
   * All operations in this class are local (no MPI).
   */
   class BondPotential 
   {

   public:

      /**
      * Constructor.
      */
      BondPotential(System& system);

      /**
      * Default constructor.
      *
      * Used only for unit testing.
      */
      BondPotential();

      /**
      * Destructor.
      */
      ~BondPotential();

      #if 0
      /**
      * Read parameters and allocate memory.
      *
      * Use iff this object was instantiated with BondPotential(System&).
      *
      * \param in input parameter stream.
      */
      virtual void readParam(std::istream& in);
      #endif

      /**
      * Create links to associated objects.
      *
      * Use iff this was instantiated with default constructor BondPotential().
      *
      * \param potential associated BondPotential object.
      */
      void associate(Boundary& boundary, 
                     const BondInteraction& interaction,
                     GroupStorage<2>& storage);

      /**
      * Add pair forces to atom forces.
      */
      void addForces();

      /**
      * Add pair forces to atom forces, and compute energy.
      */
      void addForces(double& energy);

      /**
      * Calculate total pair potential on this processor
      */
      double energy();

   private:

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated pair interaction.
      const BondInteraction* interactionPtr_;

      // Pointer to associated GroupStorage<2> object.
      GroupStorage<2>* storagePtr_;

      /**
      * Calculate bond forces and/or pair potential energy.
      */
      double addForces(bool needForce, bool needEnergy);

   };

}
#endif
