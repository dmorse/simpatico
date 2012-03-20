#ifndef PAIR_POTENTIAL_H
#define PAIR_POTENTIAL_H

#include <util/param/ParamComposite.h>  // base class
#include <ddMd/neighbor/CellList.h>     // member
#include <ddMd/neighbor/PairList.h>     // member
#include <util/boundary/Boundary.h>     // member (typedef)
#include <util/global.h>

#include <iostream>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   class System;
   class AtomStorage;
   class Domain;
   using namespace Util;

   /**
   * An PairPotential calculates forces for a parent System.
   *
   * An PairPotential has a private CellList and PairList which it
   * uses to calculate nonbonded pair forces. 
   * 
   * All operations in this class are local (no MPI).
   */
   class PairPotential : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      PairPotential(System& system);

      /**
      * Constructor (for unit testing).
      *
      * \param boundary    associated Boundary object.
      * \param domain      associated Domain object.
      * \param storage     associated AtomStorage object.
      */
      PairPotential(Boundary& boundary, Domain& domain, AtomStorage& storage);

      /**
      * Destructor.
      */
      virtual ~PairPotential();

      /**
      * Set parameters and allocate memory.
      *
      * Required if this was created with default constructor PairPotential().
      * The associate method must be called before setParam. This method sets 
      * the skin and cutoff length parameters, and allocates the internal
      * CellList and a PairList. The CellList allocates memory for the number
      * of cells required for the specified upper and lower coordinate bounds, 
      * with the specified cutoff. 
      *
      * \param lower        lower dimension of the largest expected domain.
      * \param upper        upper dimension of the largest expected domain.
      * \param cutoff       pair list skin length.
      * \param pairCapacity maximum number of pairs per processor.
      */
      void setParam(const Vector& lower, const Vector& upper,
                    double skin, int pairCapacity);

      /// \name Interaction interface
      //@{




      /**
      * Set the maximum number of atom types.
      */
      virtual void setNAtomType(int nAtomType) = 0;
  
      /**
      * Return pair energy for a single pair.
      * 
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return energy of pair
      */
      virtual 
      double energy(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Return force / separation for a single pair.
      *
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return repulsive force (< 0 if attractive) over distance
      */
      virtual 
      double forceOverR(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Return maximum cutoff.
      */
      virtual double maxPairCutoff() const = 0;

      /**
      * Return pair interaction class name (e.g., "LJPair").
      */
      virtual std::string interactionClassName() const = 0;

      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add pair forces to atom forces.
      */
      virtual void addForces() = 0;

      /**
      * Add pair forces to atom forces, and compute energy.
      */
      virtual void addForces(double& energy) = 0;

      /**
      * Calculate total pair potential on this processor
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator) = 0;
      #else
      virtual void computeEnergy();
      #endif

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
      virtual void computeStress(double& stress) const = 0;

      /**
      * Compute x, y, z bond pressure components.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const = 0;

      /**
      * Compute bond stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const = 0;
      #endif

      //@}
      /// \name Pair and Cell Lists.
      //@{

      /**
      * Build the cell and pair lists. 
      *
      * Use with objects created with PairPotential(System&). Makes a
      * CellList grid for the region defined by the associated Domain. 
      */
      void findNeighbors();

      /**
      * Build the cell and pair lists. 
      *
      * \param lower Vector of lower coordinate bounds for this processor.
      * \param lower Vector of upper coordinate bounds for this processor.
      */
      void findNeighbors(const Vector& lower, const Vector& upper);

      /**
      * Get the CellList by const reference.
      */
      const CellList& cellList() const;

      /**
      * Get the PairList by const reference.
      */
      const PairList& pairList() const;

      //@}

      /**
      * Get value of the pair list skin.
      */
      double skin() const;

      /**
      * Get value of the pair list cutoff.
      */
      double cutoff() const;

   protected:

      // CellList to construct PairList or calculate nonbonded pair forces.
      CellList cellList_;

      // Verlet pair list, to calculate nonbonded pair forces.
      PairList pairList_;

      // Boundary used to allocate space for the cell list.
      Boundary maxBoundary_;

      // Difference between pairlist cutoff and pair potential cutoff. 
      double skin_;

      // Minimum cell size = pair potential cutoff + skin.
      double cutoff_;

      // Maximum number of nonbonded pairs in pair list. 
      int pairCapacity_;

      /**
      * Read parameters and allocate memory for PairList.
      *
      * Use iff this object was instantiated with PairPotential(System&).
      *
      * \param in input parameter stream.
      */
      void readPairListParam(std::istream& in);

      /**
      * Get the PairList by const reference.
      */
      Boundary& boundary();

      /**
      * Get the PairList by const reference.
      */
      Domain& domain();

      /**
      * Get the AtomStorage by reference.
      */
      AtomStorage& storage();

   private:

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated Domain object.
      Domain* domainPtr_;

      // Pointer to associated AtomStorage object.
      AtomStorage* storagePtr_;

   };

   inline const CellList& PairPotential::cellList() const
   {  return cellList_; }

   inline const PairList& PairPotential::pairList() const
   {  return pairList_; }

   inline double PairPotential::skin() const
   {  return skin_; }

   inline double PairPotential::cutoff() const
   {  return cutoff_; }

   inline Boundary& PairPotential::boundary() 
   {  return *boundaryPtr_; }

   inline Domain& PairPotential::domain()
   {  return *domainPtr_; }

   inline AtomStorage& PairPotential::storage()
   {  return *storagePtr_; }

}
#endif
