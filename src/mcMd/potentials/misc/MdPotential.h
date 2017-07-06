#ifndef MCMD_MD_POTENTIAL_H
#define MCMD_MD_POTENTIAL_H

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
   * Potential for an MD simulation.
   *
   * \ingroup McMd_Potential_Module
   */
   class MdPotential : public ParamComposite, 
                       public EnergyCalculator, public StressCalculator
   {

   public:

      /**
      * Destructor (does nothing)
      */
      virtual ~MdPotential();

      /**
      * Add forces from this potential to all atomic forces.
      */
      virtual void addForces() = 0;

      /**
      * Return true iff this potential generates a stress contribution.
      */
      bool hasStress() const
      {  return hasStress_; }

   protected:

      /**
      * Constructor.
      *
      * Derived class constructor must set hasStress true or false.
      */
      MdPotential(bool hasStress = true);

   private:

      /// True iff potential gemerate a stress contribution.
      bool hasStress_;

   };

} 
#endif
