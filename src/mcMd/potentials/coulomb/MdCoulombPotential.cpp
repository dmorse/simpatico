#ifndef MD_COULOMB_POTENTIAL_CPP
#define MD_COULOMB_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdCoulombPotential.h" 
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MdCoulombPotential::MdCoulombPotential()
    : isInitialized_(false)
   {  setClassName("CoulombPotential"); }

   /*
   * Destructor (does nothing)
   */
   MdCoulombPotential::~MdCoulombPotential()
   {}

   void MdCoulombPotential::unsetEnergy()
   {  kSpaceEnergy_.unset(); }

   void MdCoulombPotential::unsetStress()
   {  kSpaceStress_.unset(); }

   double MdCoulombPotential::kSpaceEnergy()
   {
      if (!kSpaceEnergy_.isSet()) {
         computeEnergy();
      }
      return kSpaceEnergy_.value();
   }

   double MdCoulombPotential::rSpaceEnergy()
   {  return rSpaceAccumulator_.rSpaceEnergy(); }

   double MdCoulombPotential::energy()
   {
      double temp;
      temp  = kSpaceEnergy();
      temp += rSpaceAccumulator_.rSpaceEnergy();
      return temp;
   }

   Tensor MdCoulombPotential::kSpaceStress() 
   {
      if (!kSpaceStress_.isSet()) {
         computeStress();
      }
      return kSpaceStress_.value();
   }

   Tensor MdCoulombPotential::rSpaceStress()
   {  return rSpaceAccumulator_.rSpaceStress(); }

   Tensor MdCoulombPotential::stress()
   {
      Tensor temp;
      temp = kSpaceStress();
      temp += rSpaceAccumulator_.rSpaceStress();
      return temp;
   }
 
} 
#endif
