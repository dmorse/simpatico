#ifndef MCMD_EWALD_COULOMB_PAIR_CPP
#define MCMD_EWALD_COULOMB_PAIR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "EwaldCoulombPair.h"

namespace McMd
{

   /*
   * Set AtomTypes.
   */
   void EwaldCoulombPair::setAtomTypes(const Array<AtomType> atomTypes)
   {  atomTypesPtr_ = &atomTypes; }

}
#endif
