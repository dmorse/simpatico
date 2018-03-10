#ifndef MCMD_SPECIAL_POTENTIAL_FACADE_H
#define MCMD_SPECIAL_POTENTIAL_FACADE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/special/SpecialPotential.h> // base class

namespace McMd
{

   class System;

   using namespace Util;

   /**
   * SpecialPotential interface to an existing potential.
   *
   * \ingroup McMd_Special_Module
   */
   template <class PotentialType, class FactoryType>
   class SpecialPotentialFacade : public SpecialPotential
   {

   public:

      /**
      * Constructor.
      */
      SpecialPotentialFacade(System& system);

      /**
      * Destructor.
      */
      virtual ~SpecialPotentialFacade();

      /**
      * Read parameters.
      */
      void readParameters(std::istream& in);

      /**
      * Compute total energy.
      */
      void computeEnergy();

      /**
      * Add forces from this potential to all atomic forces.
      */
      void addForces();

   private:

      /// Pointer to parent System.
      System* systemPtr_;

      /// Pointer to potential object.
      PotentialType* potentialPtr_;

      /// Type string for potential (used by factory).
      std::string style_;

   };

}
#include "SpecialPotentialFacade.tpp" 
#endif
