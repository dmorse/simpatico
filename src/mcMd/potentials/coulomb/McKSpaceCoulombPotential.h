#ifndef MCMD_MC_KSPACE_COULOMB_POTENTIAL_H
#define MCMD_MC_KSPACE_COULOMB_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McKSpaceCoulombPotential.h"

namespace McMd
{

   class Atom;

   using namespace Util;

   /**
   * Long-range part of Coulomb potential for MC.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class McKSpaceCoulombPotential : public KSpaceCoulombPotential
   {

   public:

      /**
      * Constructor .
      */
      McKSpaceCoulombPotential();

      /**
      * Destructor (does nothing)
      */
      virtual ~McKSpaceCoulombPotential();

      /**
      * Calculate the covalent bond energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return bond potential energy of atom
      */
      virtual double atomEnergy(const Atom& atom) const = 0;

   };

} 
#endif
