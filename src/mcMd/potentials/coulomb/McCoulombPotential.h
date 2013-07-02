#ifndef MCMD_MC_COULOMB_POTENTIAL_H
#define MCMD_MC_COULOMB_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McCoulombPotential.h"

namespace McMd
{

   class Atom;

   using namespace Util;

   /**
   * Long-range part of Coulomb potential for MC.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class McCoulombPotential : public CoulombPotential
   {

   public:

      /**
      * Constructor .
      */
      McCoulombPotential();

      /**
      * Destructor (does nothing)
      */
      virtual ~McCoulombPotential();

      /**
      * Calculate the covalent bond energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return bond potential energy of atom
      */
      virtual double kspaceAtomEnergy(const Atom& atom) const = 0;

   };

} 
#endif
