#ifndef DDMD_TWO_STEP_INTEGRATOR_H
#define DDMD_TWO_STEP_INTEGRATOR_H

#include "Integrator.h"

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   class Simulation;
   using namespace Util;

   /**
   * A two-step velocity-Verlet style integrator.
   *
   * \ingroup DdMd_Integrator_Module
   */
   class TwoStepIntegrator : public Integrator
   {

   public:

      /**
      * Constructor.
      */
      TwoStepIntegrator(Simulation& simulation);

      /**
      * Destructor.
      */
      ~TwoStepIntegrator();

      /**
      * Run a simulation.
      *
      * \param nStep number of steps
      */
      void run(int nStep);

   protected:

      /**
      * Execute first step of two-step integrator.
      */
      virtual void integrateStep1() = 0;

      /**
      * Execute secodn step of two-step integrator.
      */
      virtual void integrateStep2() = 0;

   };

}
#endif
