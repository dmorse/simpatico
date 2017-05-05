#ifndef MCMD_ENERGY_CALCULATOR_H
#define MCMD_ENERGY_CALCULATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/misc/Setable.h>

namespace McMd
{

   using namespace Util;

   /**
   * Interface for a class that calculates a total energy.
   *
   * \ingroup McMd_Potential_Module
   */
   class EnergyCalculator
   {

   public:

      /**
      * Mark the energy as unknown.
      */
      virtual void unsetEnergy();

      /**
      * Calculate the total nonBonded pair energy for the associated System.
      */
      virtual void computeEnergy() = 0;

      /**
      * Return the energy contribution, compute if necessary.
      *
      * If necessary, this function calls computeEnergy() and stores the
      * value before returning stored value. 
      */
      double energy();

      //@}

   protected:

      // Setable value of energy contribution for a system.
      Setable<double> energy_;

   };

}
#endif
