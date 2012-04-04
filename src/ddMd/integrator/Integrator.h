#ifndef DDMD_INTEGRATOR_H
#define DDMD_INTEGRATOR_H

#include <util/param/ParamComposite.h>          // base class
#include <ddMd/simulation/SimulationAccess.h>   // base class
#include <ddMd/util/DdTimer.h>                  // member

#include <iostream>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * An Integrator numerically integrates the equations of motion. 
   *
   * \ingroup DdMd_Integrator_Module
   */
   class Integrator : public ParamComposite, public SimulationAccess
   {

   public:

      /**
      * Constructor.
      */
      Integrator(Simulation& simulation);

      /**
      * Destructor.
      */
      ~Integrator();

      /**
      * Initialize just before integration.
      */
      virtual void setup() = 0;

      /**
      * Implement one step.
      */
      virtual void step() = 0;

      /**
      * Run a simulation of iStep steps.
      */
      virtual void run(int nStep) = 0;

      /**
      * Output statistics immediately after a run.
      */
      virtual void outputStatistics(std::ostream& out);

   protected:

      enum TimeId {DIAGNOSTIC, INTEGRATE1, EXCHANGE, NEIGHBOR, UPDATE, 
                   FORCE, INTEGRATE2, NTime};

      DdTimer& timer()
      { return timer_; }

      int   nStep_;

      int   iStep_;

   private:

      DdTimer timer_;

   };


}
#endif
