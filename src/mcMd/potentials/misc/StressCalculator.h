#ifndef MCMD_STRESS_CALCULATOR_H
#define MCMD_STRESS_CALCULATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Tensor.h>
#include <util/misc/Setable.h>

namespace Util {
   class Vector;
}

namespace McMd
{

   using namespace Util;

   /**
   * Interface for a stress calculator.
   *
   * \ingroup McMd_Potential_Module
   */
   class StressCalculator
   {

   public:

      /**
      * Mark the stress as unknown.
      */
      virtual void unsetStress();

      /**
      * Compute and store the stress tensor.
      *
      * Default implementation throws an Exception, to allow 
      * testing and graceful failure potentials that do not 
      * create stress.
      */
      virtual void computeStress()
      {  UTIL_THROW("Unimplemented computeStress function"); }
      

      /**
      * Get pair stress tensor.
      *
      * If necessary, this function calls computeStress() before
      * accessing value.
      *
      * \param stress (output) pair stress tensor
      */
      void computeStress(Tensor& stress);

      /**
      * Get the xx, yy, zz non-Coulomb pair pressures.
      *
      * If necessary, this function calls computeStress() before
      * accessing values.
      *
      * \param pressures (output) diagonal pair stress components
      */
      void computeStress(Vector& pressures);

      /**
      * Get the scalar pressure.
      *
      * If necessary, this function calls computeStress() before
      * accessing values.
      *
      * \param pressure (output) scalar pair pressure.
      */
      void computeStress(double& pressure);

      /**
      * Return false if subclass does not generate stress.
      */
      bool createsStress() const;

   protected:

      // Setable value of stress tensor.
      Setable<Tensor> stress_;

      /**
      * Constructor (protected to prevent direct instantiation).
      *
      * Derived class constructor must set hasStress true or false.
      */
      StressCalculator(bool createsStress = true);

   private:

      /// True iff potential gemerate a stress contribution.
      bool createsStress_;

   };

}
#endif
