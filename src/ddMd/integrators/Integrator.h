#ifndef DDMD_INTEGRATOR_H
#define DDMD_INTEGRATOR_H

#include <util/param/ParamComposite.h>          // base class
#include <ddMd/simulation/SimulationAccess.h>   // base class
#include <ddMd/util/DdTimer.h>                  // member

#include <iostream>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
      * Run a simulation of iStep steps.
      */
      virtual void run(int nStep) = 0;

      /**
      * Output statistics immediately after a run.
      */
      virtual void outputStatistics(std::ostream& out);

      /**
      * Get time per processor of previous run.  
      */
      double time() const;

      /**
      * Get number of steps in last run.
      */ 
      int nStep() const;

   protected:

      /// Timestamps for loop timing.
      enum TimeId {DIAGNOSTIC, INTEGRATE1, CHECK, TRANSFORM_F, EXCHANGE, 
                   CELLLIST, TRANSFORM_R, PAIRLIST, UPDATE, PAIR_FORCE, 
                   BOND_FORCE, ANGLE_FORCE, DIHEDRAL_FORCE,
                   INTEGRATE2, NTime};

      /// Return the timer by reference.
      DdTimer& timer()
      {  return timer_; }

      /**
      * Setup state of atoms just before integration.
      * 
      * Should be called in all subclass setup methods.
      */
      void setupAtoms();

      /**
      * Compute forces for all local atoms, with timing.
      *
      * Identical to Simulation::forceCompute(), with added timing.
      * Upon return, forces are correct for all local atoms. Values
      * of the forces on ghost atoms are undefined. Executes reverse
      * communication if needed, and emits Simulation::forceSignal().
      */
      void computeForces();

      /// Total number of steps
      int   nStep_;

      /// Current step number.
      int   iStep_;

   private:

      // Performance timer
      DdTimer timer_;

   };

}
#endif
