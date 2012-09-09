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
#include <util/util/Log.h>
#include <util/global.h>

// Uncomment to enable paranoid validity checks.
//#define  DDMD_INTEGRATOR_DEBUG

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
      bool needExchange;
      nStep_ = nStep;


      // Compute nAtomTotal.
      atomStorage().unsetNAtomTotal();
      atomStorage().computeNAtomTotal(domain().communicator());

      // Clear all stored computations, compute nAtomTotal.
      simulation().modifySignal().notify();
      clearStatistics();

      // Main MD loop
      timer().start();
      exchanger().timer().start();
      for (iStep_ = 0; iStep_ < nStep_; ++iStep_) {

         if (Diagnostic::baseInterval > 0) {
            if (iStep_ % Diagnostic::baseInterval == 0) {
               simulation().diagnosticManager().sample(iStep_);
            }
         }
         timer().stamp(DIAGNOSTIC);
 
         // First step of integration: Update positions, half velocity 
         integrateStep1();
         timer().stamp(INTEGRATE1);
  
         // Unset precomputed values of all energies, stresses, etc.
         simulation().modifySignal().notify();

         #ifdef DDMD_INTEGRATOR_DEBUG
         simulation().isValid();
         #endif

         // Check if exchange and reneighboring is necessary
         needExchange = atomStorage().needExchange(domain().communicator(), 
                                                   pairPotential().skin());
         timer().stamp(Integrator::CHECK);

         // Exchange atoms if necessary
         if (needExchange) {
            atomStorage().clearSnapshot();
            if (!UTIL_ORTHOGONAL && atomStorage().isCartesian()) {
               atomStorage().transformCartToGen(boundary());
               timer().stamp(Integrator::TRANSFORM_F);
            }
            exchanger().exchange();
            timer().stamp(Integrator::EXCHANGE);
            pairPotential().buildCellList();
            timer().stamp(Integrator::CELLLIST);
            if (!UTIL_ORTHOGONAL) {
               atomStorage().transformGenToCart(boundary());
               timer().stamp(Integrator::TRANSFORM_R);
            }
            atomStorage().makeSnapshot();
            pairPotential().buildPairList();
            timer().stamp(Integrator::PAIRLIST);
         } else {
            exchanger().update();
            timer().stamp(UPDATE);
         }
   
         #ifdef DDMD_INTEGRATOR_DEBUG
         simulation().isValid();
         #endif

         // Calculate new forces for all local atoms. Also calculate virial stresses
         // for constant pressure ensemble. Both methods send the modifyForce signal.
         if (simulation().boundaryEnsemble().isRigid()) {
            computeForces();
         } else {
            computeForcesAndVirial();
         }

         // 2nd step of integration, which finishes the velocity update.
         // This method normally call simulation().velocitySignal().notify()
         integrateStep2();
         timer().stamp(INTEGRATE2);
   
         #ifdef DDMD_INTEGRATOR_DEBUG
         simulation().isValid();
         #endif

      }
      exchanger().timer().stop();
      timer().stop();

      // Final diagnostic output
      simulation().diagnosticManager().output();
             
      // Reduce statistics for run
      #ifdef UTIL_MPI
      timer().reduce(domain().communicator());
      exchanger().timer().reduce(domain().communicator());
      #endif

      if (domain().isMaster()) {
         Log::file() << std::endl;
      }

   }

}
#endif
