#ifndef DDMD_TWO_STEP_INTEGRATOR_CPP
#define DDMD_TWO_STEP_INTEGRATOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TwoStepIntegrator.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/potentials/pair/PairPotential.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <util/misc/Log.h>
#include <util/global.h>

// Uncomment to enable paranoid validity checks.
//#define DDMD_INTEGRATOR_DEBUG

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   TwoStepIntegrator::TwoStepIntegrator(Simulation& simulation)
    : Integrator(simulation)
   {}

   /*
   * Destructor.
   */
   TwoStepIntegrator::~TwoStepIntegrator()
   {}

   /*
   * Run integrator for nStep steps.
   */
   void TwoStepIntegrator::run(int nStep)
   {
      // Precondition
      if (atomStorage().isCartesian()) {
         UTIL_THROW("Error: Atom coordinates are Cartesian");
      }

      // Unset all stored computations.
      simulation().modifySignal().notify();

      // Recompute nAtomTotal.
      atomStorage().unsetNAtomTotal();
      atomStorage().computeNAtomTotal(domain().communicator());

      // Setup required before main loop (ghosts, forces, etc)
      // Atomic coordinates are Cartesian on exit from setup.
      setup();
      simulation().diagnosticManager().setup();

      // Main MD loop
      timer().start();
      exchanger().timer().start();
      int  beginStep = iStep_;
      int  endStep = iStep_ + nStep;
      bool needExchange;
      for ( ; iStep_ < endStep; ++iStep_) {

         // Atomic coordinates must be Cartesian on entry to loop body.
         if (!atomStorage().isCartesian()) {
            UTIL_THROW("Error: Atomic coordinates are not Cartesian");
         }

         // Sample diagnostics, if scheduled.
         simulation().diagnosticManager().sample(iStep_);

         // Write restart file, if scheduled.
         if (saveInterval() > 0) {
            if (iStep_ % saveInterval() == 0) {
               if (iStep_ > beginStep) {
                  simulation().save(saveFileName());
               }
            }
         }
         timer().stamp(DIAGNOSTIC);
 
         // First step of integration: Update positions, half velocity 
         integrateStep1();
         timer().stamp(INTEGRATE1);
  
         // Unset precomputed values of all energies, stresses, etc.
         simulation().modifySignal().notify();

         #ifdef DDMD_INTEGRATOR_DEBUG
         // Sanity check
         simulation().isValid();
         #endif

         // Check if exchange and reneighboring is necessary
         needExchange = isExchangeNeeded(pairPotential().skin());

         if (!atomStorage().isCartesian()) {
            UTIL_THROW("Error: atomic coordinates are not Cartesian");
         }

         // Exchange atoms if necessary
         if (needExchange) {
            atomStorage().clearSnapshot();
            atomStorage().transformCartToGen(boundary());
            timer().stamp(Integrator::TRANSFORM_F);
            exchanger().exchange();
            timer().stamp(Integrator::EXCHANGE);
            pairPotential().buildCellList();
            timer().stamp(Integrator::CELLLIST);
            atomStorage().transformGenToCart(boundary());
            timer().stamp(Integrator::TRANSFORM_R);
            atomStorage().makeSnapshot();
            pairPotential().buildPairList();
            timer().stamp(Integrator::PAIRLIST);
         } else {
            exchanger().update();
            timer().stamp(UPDATE);
         }
         simulation().exchangeSignal().notify();
   
         #ifdef DDMD_INTEGRATOR_DEBUG
         // Sanity check
         simulation().isValid();
         #endif

         if (!atomStorage().isCartesian()) {
            UTIL_THROW("Error: Atomic coordinates are not Cartesian");
         }

         // Calculate new forces for all local atoms. If constant pressure
         // ensembles (not rigid), calculate virial stress. Both methods 
         // send the modifyForce signal.
         if (simulation().boundaryEnsemble().isRigid()) {
            computeForces();
         } else {
            computeForcesAndVirial();
         }

         // 2nd step of integration. This finishes the velocity update.
         // This method normally calls simulation().velocitySignal().notify()
         integrateStep2();
         timer().stamp(INTEGRATE2);
   
         #ifdef DDMD_INTEGRATOR_DEBUG
         // Sanity check
         simulation().isValid();
         #endif

      }
      exchanger().timer().stop();
      timer().stop();

      // Final diagnostics and restart file, if scheduled.
      simulation().diagnosticManager().sample(iStep_);
      if (saveInterval() > 0) {
         if (iStep_ % saveInterval() == 0) {
            simulation().save(saveFileName());
         }
      }

      // Transform to scaled coordinates, in preparation for the next run.
      atomStorage().transformCartToGen(boundary());

      if (domain().isMaster()) {
         Log::file() << std::endl;
      }

   }

}
#endif
