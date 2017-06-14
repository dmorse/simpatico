#ifndef MCMD_MD_PAIR_POTENTIAL_H
#define MCMD_MD_PAIR_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/pair/PairPotential.h>   // base class
#include <mcMd/simulation/SystemInterface.h>      // base class
#include <mcMd/neighbor/PairList.h>               // member

#include <util/param/ParamComposite.h>            // base class
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
   * \ingroup McMd_Pair_Module
   */
   class MdPairPotential : public ParamComposite, public PairPotential,
                           public SystemInterface
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

      /**
      * Calculate non-bonded pair forces for all atoms in this System.
      *
      * Adds non-bonded pair forces to the current values of the forces
      * for all atoms in this system. Before calculating forces, this
      * function checks if the pair list is current, and rebuilds it if 
      * necessary.
      */
      virtual void addForces() = 0;

      /// \name PairList Management (non-virtual)
      //@{

      /**
      * Build the internal PairList.
      *
      * Before building the PairList, this method also shifts all atom
      * positions into the primary image of the periodic boundary. In 
      * an MD simulation, the positions are allowed to wander out of 
      * the primary image between calls to this function.
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

   };

   // Inline functions
  
   /*
   * Return PairList by reference.
   */
   inline const PairList& MdPairPotential::pairList() const
   {  return pairList_; }

} 
#endif
