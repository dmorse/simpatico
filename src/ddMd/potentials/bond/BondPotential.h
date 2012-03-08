#ifndef BOND_POTENTIAL_H
#define BOND_POTENTIAL_H

#include <util/param/ParamComposite.h>        // base class
#include <ddMd/boundary/Boundary.h>           // typedef

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
   class BondPotential  : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      BondPotential(System& system);

      /**
      * Constructor (for unit testing).
      *
      * \param boundary associated Boundary object.
      * \param storage  associated bond storage.
      */
      BondPotential(Boundary& boundary, GroupStorage<2>& storage);

      /**
      * Destructor.
      */
      ~BondPotential();

      /// \name Interaction interface
      //@{

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNBondType(int nBondType) = 0;
  
      /**
      * Return pair energy for a single pair.
      */
      virtual double energy(double rsq, int bondTypeId) const = 0;

      /**
      * Return force / separation for a single pair.
      */
      virtual double forceOverR(double rsq, int bondTypeId) const = 0;

      /**
      * Return force / separation for a single pair.
      */
      virtual double 
      randomBondLength(Random* random, double beta, int bondTypeId) const = 0;

      #if 0
      /**
      * Return pair interaction class name (e.g., "HarmonicBond").
      */
      virtual std::string interactionClassName() const = 0;
      #endif

      //@}

      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add the bond forces for all atoms.
      */
      virtual void addForces() = 0;

      /**
      * Add pair forces to atom forces, and compute energy.
      */
      virtual void addForces(double& energy) = 0;

      /**
      * Calculate total pair potential on this processor
      */
      virtual double energy() = 0;

      #if 0
      /**
      * Compute total bond pressure.
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const;

      /**
      * Compute x, y, z bond pressure components.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const;

      /**
      * Compute bond stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const;

      //@}
      #endif

   protected:

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated GroupStorage<2> object.
      GroupStorage<2>* storagePtr_;

   };

}
#endif
