#ifndef DDMD_BOND_POTENTIAL_H
#define DDMD_BOND_POTENTIAL_H

#include <ddMd/potentials/Potential.h>        // base class
#include <util/boundary/Boundary.h>           // typedef

#include <iostream>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   using namespace Util;

   class Simulation;
   template <int N> class GroupStorage;

   /**
   * Calculates bond forces and energies for a parent Simulation.
   *
   * \ingroup DdMd_Bond_Module
   */
   class BondPotential  : public Potential
   {

   public:

      /**
      * Constructor.
      */
      BondPotential(Simulation& simulation);

      /**
      * Default constructor (for unit testing).
      */
      BondPotential();

      /**
      * Associate with related objects.
      *
      * Call iff object instantiated with default constructor.
      *
      * \param boundary associated Boundary object.
      * \param storage  associated GroupStorage<2> object.
      */
      void associate(Boundary& boundary, GroupStorage<2>& storage);

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
      virtual double bondEnergy(double rsq, int bondTypeId) const = 0;

      /**
      * Return force / separation for a single pair.
      */
      virtual double bondForceOverR(double rsq, int bondTypeId) const = 0;

      /**
      * Return force / separation for a single pair.
      */
      virtual double 
      randomBondLength(Random* random, double beta, int bondTypeId) const = 0;

      /**
      * Modify a bond interaction parameter, identified by a string.
      *
      * \param name       parameter variable name
      * \param bondTypeId bond type index
      * \param value      new value of parameter
      */
      virtual void set(std::string name, int bondTypeId, double value) = 0;

      /**
      * Get a bond parameter value, identified by a string.
      *
      * \param name       parameter variable name
      * \param bondTypeId bond type index
      */
      virtual double get(std::string name, int bondTypeId) const = 0;

      /**
      * Return pair interaction class name (e.g., "HarmonicBond").
      */
      virtual std::string interactionClassName() const = 0;

      //@}

   protected:

      /**
      *  Return boundary by reference.   
      */
      Boundary& boundary() const;

      /**
      *  Return bond storage by reference.   
      */
      GroupStorage<2>& storage() const;

   private:

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated GroupStorage<2> object.
      GroupStorage<2>* storagePtr_;

   };

   inline Boundary& BondPotential::boundary() const
   { return *boundaryPtr_; }

   inline GroupStorage<2>& BondPotential::storage() const
   { return *storagePtr_; }

}
#endif
