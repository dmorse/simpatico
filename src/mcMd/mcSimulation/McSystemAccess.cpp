/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSystemAccess.h"
#include <mcMd/mcSimulation/McSystem.h>  

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   *
   * \param system parent McSystem
   */
   McSystemAccess::McSystemAccess(McSystem& mcSystem)
    : SubSystem(mcSystem)
      #ifndef INTER_NOPAIR
      , pairPotentialPtr_(0)
      #endif
      #ifdef INTER_BOND
      , bondPotentialPtr_(0)
      #endif
      #ifdef INTER_ANGLE
      , anglePotentialPtr_(0)
      #endif
      #ifdef INTER_DIHEDRAL
      , dihedralPotentialPtr_(0)
      #endif
      #ifdef INTER_EXTERNAL
      , externalPotentialPtr_(0)
      #endif
   {
      #ifndef INTER_NOPAIR
      pairPotentialPtr_ = &mcSystem.pairPotential();
      #endif
      #ifdef INTER_BOND
      if (hasBonds()) {
         bondPotentialPtr_ = &mcSystem.bondPotential();
      }
      #endif
      #ifdef INTER_ANGLE
      if (hasAngles()) {
         anglePotentialPtr_ = &mcSystem.anglePotential();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedrals()) {
         dihedralPotentialPtr_ = &mcSystem.dihedralPotential();
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal()) {
         externalPotentialPtr_ = &mcSystem.externalPotential();
      }
      #endif
   }
   
   /*
   * Destructor.
   */
   McSystemAccess::~McSystemAccess()
   {}

}
