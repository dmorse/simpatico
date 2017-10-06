#ifndef MCMD_SPECIAL_POTENTIAL_H
#define MCMD_SPECIAL_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>               // base class
#include <mcMd/potentials/misc/EnergyCalculator.h>   // base class
#include <mcMd/potentials/misc/StressCalculator.h>   // base class

namespace McMd
{

   using namespace Util;

   /**
   * Specialized potential for an MD simulation.
   *
   * \ingroup McMd_Potential_Module
   */
   class SpecialPotential : public ParamComposite, 
                            public EnergyCalculator, public StressCalculator
   {

   public:

      /**
      * Destructor (does nothing)
      */
      virtual ~SpecialPotential();

      /**
      * Add forces from this potential to all atomic forces.
      */
      virtual void addForces() = 0;

   protected:

      /**
      * Constructor.
      *
      * Derived class constructor must set createsStress true or false.
      */
      SpecialPotential(bool createsStress = true);

   };

} 
#endif
