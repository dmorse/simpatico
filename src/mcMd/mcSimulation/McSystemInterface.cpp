/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSystemInterface.h"
#include <mcMd/mcSimulation/McSystem.h>  

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   *
   * \param system parent McSystem
   */
   McSystemInterface::McSystemInterface(McSystem& mcSystem)
    : SystemInterface(mcSystem)
      #ifndef SIMP_NOPAIR
      , pairPotentialPtr_(0)
      #endif
      #ifdef SIMP_BOND
      , bondPotentialPtr_(0)
      #endif
      #ifdef SIMP_ANGLE
      , anglePotentialPtr_(0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , dihedralPotentialPtr_(0)
      #endif
      #ifdef SIMP_EXTERNAL
      , externalPotentialPtr_(0)
      #endif
   {
      #ifndef SIMP_NOPAIR
      pairPotentialPtr_ = &mcSystem.pairPotential();
      #endif
      #ifdef SIMP_BOND
      if (hasBonds()) {
         bondPotentialPtr_ = &mcSystem.bondPotential();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (hasAngles()) {
         anglePotentialPtr_ = &mcSystem.anglePotential();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedrals()) {
         dihedralPotentialPtr_ = &mcSystem.dihedralPotential();
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal()) {
         externalPotentialPtr_ = &mcSystem.externalPotential();
      }
      #endif
   }
   
   /*
   * Destructor.
   */
   McSystemInterface::~McSystemInterface()
   {}

}
