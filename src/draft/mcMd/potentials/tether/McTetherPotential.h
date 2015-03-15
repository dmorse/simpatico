#ifdef  INTER_TETHER
#ifndef MCMD_MC_TETHER_POTENTIAL_H
#define MCMD_MC_TETHER_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/tether/TetherPotential.h>  // base class

namespace McMd
{

   using namespace Util;

   class System;

   /**
   * A TetherPotential for MC simulations.
   *
   * \ingroup McMd_Tether_Module
   */
   class McTetherPotential : public TetherPotential 
   {

   #if 0
   public:

      /**   
      * Constructor.
      */
      McTetherPotential();

      /** 
      * Destructor.
      */
      virtual ~McTetherPotential();

      /**
      * Calculate the tether energy for a specified Atom.
      *
      * \param  atom Atom object of interest
      * \return tether potential energy of atom
      */
      virtual double atomEnergy(const Atom& atom) const = 0;
   #endif

   };

} 
#endif
#endif
