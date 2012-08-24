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
#include <util/util/Log.h>
#include <util/global.h>

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
         simulation().modifySignal().notify();
         timer().stamp(INTEGRATE1);
   
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
   
         // Calculate new forces for all local atoms
         computeForces();
         simulation().forceSignal().notify();

         // 2nd step of integration: Finish velocity update.
         integrateStep2();
         simulation().velocitySignal().notify();
         timer().stamp(INTEGRATE2);
   
      }
      exchanger().timer().stop();
      timer().stop();

      // Compute and reduce statistics for run
      #ifdef UTIL_MPI
      timer().reduce(domain().communicator());
      exchanger().timer().reduce(domain().communicator());
      pairPotential().pairList().computeStatistics(domain().communicator());
      atomStorage().computeNAtomTotal(domain().communicator());
      #else
      pairPotential().pairList().computeStatistics();
      #endif

      if (domain().isMaster()) {
         Log::file() << std::endl;
      }

   }

}
#endif
