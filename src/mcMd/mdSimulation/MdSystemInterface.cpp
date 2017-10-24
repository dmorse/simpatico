/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdSystemInterface.h"
#include <mcMd/mdSimulation/MdSystem.h>  

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   *
   * \param system parent MdSystem
   */
   MdSystemInterface::MdSystemInterface(MdSystem& mdSystem)
    : SystemInterface(mdSystem)
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
      pairPotentialPtr_ = &mdSystem.pairPotential();
      #endif
      #ifdef SIMP_BOND
      if (hasBonds()) {
         bondPotentialPtr_ = &mdSystem.bondPotential();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (hasAngles()) {
         anglePotentialPtr_ = &mdSystem.anglePotential();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedrals()) {
         dihedralPotentialPtr_ = &mdSystem.dihedralPotential();
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal()) {
         externalPotentialPtr_ = &mdSystem.externalPotential();
      }
      #endif
   }
   
   /*
   * Destructor.
   */
   MdSystemInterface::~MdSystemInterface()
   {}

}
