#ifndef MCMD_NOPAIR
#ifndef MD_PAIR_POTENTIAL_H
#define MD_PAIR_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/SubSystem.h>             // base class
#include <util/param/ParamComposite.h>             // base class
#include <mcMd/potentials/pair/PairPotential.h>  // base class
#include <mcMd/neighbor/PairList.h>                // member

#include <util/global.h>

namespace Util
{
   class Vector;
   class Tensor;
}

namespace McMd
{

   using namespace Util;

   class System;

   /**
   * An PairPotential for MD simulation.
   *
   * \ingroup Pair_Module
   */
   class MdPairPotential : public ParamComposite, public SubSystem, public PairPotential
   {

   public:

      /**   
      * Constructor.
      */
      MdPairPotential(System& system);

      /** 
      * Destructor.
      */
      virtual ~MdPairPotential();

      /// \name Energy and Force Evaluators (pure virtual)
      //@{

      #if 0
      /**
      * Return pair energy for a single pair.
      */
      virtual double energy(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Return force / separation for a single pair.
      */
      virtual double forceOverR(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Return maximum cutoff.
      */
      virtual double maxPairCutoff() const = 0;

      /**
      * Calculate the total nonBonded pair energy for this System.
      * 
      * Rebuilds the PairList if necessary before calculating energy.
      */
      virtual double energy() = 0;

      /**
      * Compute total nonbonded pressure
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const = 0;

      /**
      * Compute x, y, z nonbonded pressures.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const = 0;

      /**
      * Compute stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const = 0;
      #endif

      /**
      * Calculate non-bonded pair forces for all atoms in this System.
      *
      * Adds non-bonded pair forces to the current values of the
      * forces for all atoms in this system. Before calculating 
      * forces, the method checks if the pair list is current, 
      * and rebuilds it if necessary.
      */
      virtual void addForces() = 0;

      //@}
      /// \name PairList Management (non-virtual)
      //@{

      /**
      * Build the internal PairList.
      *
      * As part of building the PairList, this method also shifts all atom
      * positions into the primary image of the periodic boundary. In an
      * MD simulation, the positions are allowed to wander out of the
      * primary image between calls to this function.
      */
      void buildPairList();

      /**
      * Return true if PairList is current, false if obsolete.
      *
      * A PairList is considered current iff it has been built
      * previously, and no particle has moved by a distance 
      * greater than skin/2 since the PairList was last built.
      */
      bool isPairListCurrent();

      /**
      * Clear all statistical accumulators stored in PairList.
      */
      void clearPairListStatistics();

      /**
      * Return a const reference to the internal PairList
      */
      const PairList& pairList() const;

      //@}

   protected:

      /// Verlet neighbor pair list for nonbonded interactions.
      PairList pairList_;

      /// Maximum boundary (used to allocated memory for PairList).
      Boundary maxBoundary_;

   };

   // Inline functions
  
   /*
   * Return PairList by reference.
   */
   inline const PairList& MdPairPotential::pairList() const
   {  return pairList_; }

} 
#endif
#endif
