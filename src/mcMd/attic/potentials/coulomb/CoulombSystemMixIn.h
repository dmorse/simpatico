#ifndef MCMD_COULOMB_SYSTEM_MIXIN_H
#define MCMD_COULOMB_SYSTEM_MIXIN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/chemistry/AtomType.h>   // Member template parameter
#include <util/containers/Array.h>     // member class template
#include <util/boundary/Boundary.h>    // typedef

namespace McMd
{

   class Simulation;
   class System;

   using namespace Util;

   /**
   * MixIn class for implementation of CoulombPotential subclasses.
   *
   * Concrete subclasses of CoulombPotential should also be derived from
   * CoulombSystemMixIn, as a private or protected base class. This class
   * simply provides pointers to the parent System, Simulation, Boundary, 
   * and an array of AtomType objects, for use in implementation of the
   * kspace computation member functions.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class CoulombSystemMixIn
   {

   public:

      /**
      * Constructor.
      */
      CoulombSystemMixIn(System& system);

   protected:

      /// Pointer to associated Simulation
      Simulation* simulationPtr_;

      /// Pointer to associated System
      System* systemPtr_;

      /// Pointer to associated Boundary
      Boundary* boundaryPtr_;

      /// Pointer to array of AtomTypes.
      const Array<AtomType>* atomTypesPtr_;

   };

} 
#endif
