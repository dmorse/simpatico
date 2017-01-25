#ifndef DDMD_BOND_POTENTIAL_H
#define DDMD_BOND_POTENTIAL_H

#include <ddMd/potentials/Potential.h>        // base class
#include <util/boundary/Boundary.h>           // typedef

#include <iostream>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   using namespace Util;

   class Simulation;
   template <int N> class GroupStorage;

   /**
   * Abstract base class for computing bond forces and energies.
   *
   * \ingroup DdMd_Bond_Module
   */
   class BondPotential  : public Potential
   {

   public:

      /**
      * Constructor.
      *
      * This is the constructor that is used during a simulation.
      *
      * \param simulation  parent Simulation object
      */
      BondPotential(Simulation& simulation);

      /**
      * Default constructor.
      *
      * This constructor is provided only to simplify unit testing.
      */
      BondPotential();

      /**
      * Destructor.
      */
      ~BondPotential();

      /**
      * Create association with related objects.
      *
      * Call iff object instantiated with default constructor, for
      * unit testing.
      *
      * \param boundary  associated Boundary object.
      * \param storage  associated GroupStorage<2> object.
      */
      void associate(Boundary& boundary, GroupStorage<2>& storage);

      /**
      * Set the maximum number of atom types.
      *
      * The implementation by a subclass should set the nBondType of the
      * associated Interaction. This should be called by the main constructor
      * of a concrete subclass, to set nBondType to simulation.nBondType(), 
      * or by a user if the object is instantiated with default constructor 
      * for unit testing.
      *
      * \param nBondType  maximum number of bond types (max index + 1).
      */
      virtual void setNBondType(int nBondType) = 0;
  
      /// \name Interaction interface
      //@{

      /**
      * Compute and return pair energy for a single pair.
      *
      * \param rsq  square of distance between atoms
      * \param bondTypeId  bond type index
      */
      virtual double bondEnergy(double rsq, int bondTypeId) const = 0;

      /**
      * Compute and return force / separation for a single pair.
      *
      * \param rsq  square of distance between atoms
      * \param bondTypeId  bond type index
      */
      virtual double bondForceOverR(double rsq, int bondTypeId) const = 0;

      /**
      * Return a random bond length, chosen from a Boltzmann distribution.
      *
      * This function should return the length of a random bond vector 
      * from the Boltzmann distribution for the specified bond interaction,
      * giving a probability distribution r^2 exp( -U(r) / kT), where U(r)
      * is the bond potential.  
      *
      * \param random  pointer to a random number generator
      * \param beta  inverse thermal energy 1/T, with T in units of energy
      * \param bondTypeId  bond type index
      */
      virtual double 
      randomBondLength(Random* random, double beta, int bondTypeId) const = 0;

      /**
      * Modify a bond interaction parameter, identified by a string.
      *
      * \param name  parameter variable name
      * \param bondTypeId  bond type index
      * \param value  new value of parameter
      */
      virtual void set(std::string name, int bondTypeId, double value) = 0;

      /**
      * Get a bond parameter value, identified by a string.
      *
      * \param name  interaction parameter name
      * \param bondTypeId  bond type index
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
