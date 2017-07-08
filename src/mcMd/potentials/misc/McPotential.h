#ifndef MCMC_MC_POTENTIAL_H
#define MCMC_MC_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>               // base class
#include <mcMd/potentials/misc/EnergyCalculator.h>   // base class
#include <mcMd/potentials/misc/StressCalculator.h>   // base class

#include <string>

namespace Util
{
   class Vector;
   class Tensor;
   class Random;
}

namespace McMd
{

   using namespace Util;

   class Atom;

   /**
   * Potential for an MC simulation.
   *
   * \ingroup McMd_Potential_Module
   */
   class McPotential : public virtual ParamComposite, 
                       public virtual EnergyCalculator, 
                       public virtual StressCalculator
   {

   public:

      /**
      * Destructor (does nothing)
      */
      virtual ~McPotential();

      /**
      * Compute and return energy for one atom.
      */
      virtual double atomEnergy(const Atom& atom) const = 0;

   protected:

      /**
      * Constructor.
      *
      * Derived class constructor must set hasStress true or false.
      */
      McPotential(bool createsStress = true);

   };

} 
#endif
