#ifndef PAIR_POTENTIAL_H
#define PAIR_POTENTIAL_H

#include <util/param/ParamComposite.h>
#include <ddMd/potentials/PairInteraction.h>
#include <ddMd/neighbor/CellList.h>
#include <ddMd/neighbor/PairList.h>

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
      * Default constructor.
      *
      * Used only for unit testing.
      */
      PairPotential();

      /**
      * Destructor.
      */
      ~PairPotential();

      /**
      * Read parameters and allocate memory.
      *
      * Use iff this object was instantiated with PairPotential(System&).
      *
      * \param in input parameter stream.
      */
      virtual void readParam(std::istream& in);

      /**
      * Create links to associated objects.
      *
      * Use iff this was instantiated with default constructor PairPotential().
      *
      * \param storage   associated AtomStorage object.
      * \param potential associated PairInteraction object.
      */
      void associate(AtomStorage& storage, const PairInteraction& potential);

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

      /**
      * Set integer id to specify algorithm.
      *
      * \param methodId algorithm id: 0=pair list, 1=cell list, 2=N^2 loop.
      */
      void setMethodId(int methodId);
  
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
      * Set forces for all local atoms to zero.
      */
      void zeroForces();

      /**
      * Calculate forces for all local atoms.
      */
      void calculateForces();

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

      /**
      * Get value of the pair list skin.
      */
      double skin() const;

      /**
      * Get value of the pair list cutoff.
      */
      double cutoff() const;

      /**
      * Get the CellList by const reference.
      */
      const CellList& cellList() const;

      /**
      * Get the PairList by const reference.
      */
      const PairList& pairList() const;

   private:

      // CellList to construct PairList or calculate nonbonded pair forces.
      CellList cellList_;

      // Verlet pair list, to calculate nonbonded pair forces.
      PairList pairList_;

      // Boundary used to allocate space for the cell list.
      Boundary maxBoundary_;

      // Difference between pairlist cutoff and pair potential cutoff. 
      double   skin_;

      // Minimum cell size = pair potential cutoff + skin.
      double   cutoff_;

      // Pointer to associated AtomStorage object.
      AtomStorage* storagePtr_;

      // Pointer to associated pair interaction.
      const PairInteraction* interactionPtr_;

      // Pointer to associated Domain object.
      const Domain* domainPtr_;

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Maximum number of nonbonded pairs in pair list. 
      int pairCapacity_;

      // Index for method used to calculate forces / energies.
      int methodId_;

      /**
      * Calculate atomic pair forces and/or pair potential energy.
      *
      * Use the Verlet pair list. 
      */
      double addForcesList(bool needForce, bool needEnergy);

      /**
      * Calculate atomic pair forces and/or pair potential energy.
      *
      * Use cell list (but not pair list).
      */
      double addForcesCell(bool needForce, bool needEnergy);

      /**
      * Calculate atomic pair forces and/or pair potential energy.
      * 
      * Use an O(N^2) double loop over all atoms.
      */
      double addForcesNSq(bool needForce, bool needEnergy);

   };

   inline const CellList& PairPotential::cellList() const
   {  return cellList_; }

   inline const PairList& PairPotential::pairList() const
   {  return pairList_; }

   inline double PairPotential::skin() const
   {  return skin_; }

   inline double PairPotential::cutoff() const
   {  return cutoff_; }

}
#endif
