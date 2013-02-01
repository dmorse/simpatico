#ifndef MCMD_DIHEDRAL_POTENTIAL_CPP
#define MCMD_DIHEDRAL_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DihedralPotential.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DihedralPotential::DihedralPotential()
   {  setClassName("DihedralPotential"); }

   /*
   * Destructor (does nothing)
   */
   DihedralPotential::~DihedralPotential()
   {}

} 
#endif
