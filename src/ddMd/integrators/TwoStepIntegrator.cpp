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
#ifdef DDMD_MODIFIERS
#include <ddMd/modifiers/ModifierManager.h>
#endif
#include <ddMd/analyzers/AnalyzerManager.h>
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

      // References to managers
      #ifdef DDMD_MODIFIERS 
      ModifierManager& modifierManager = simulation().modifierManager();
      #endif
      AnalyzerManager& analyzerManager = simulation().analyzerManager();

      // Unset all stored computations.
      simulation().modifySignal().notify();

      // Recompute nAtomTotal.
      atomStorage().unsetNAtomTotal();
      atomStorage().computeNAtomTotal(domain().communicator());

      // Setup required before main loop (atoms, ghosts, groups, forces, etc)
      // Atomic coordinates are Cartesian on exit from Integrator::setup().
      // Atomic forces should be set on exit from Integrator::setup().
      setup();
      #ifdef DDMD_MODIFIERS 
      modifierManager.setup();
      #endif
      analyzerManager.setup();

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

         // Sample analyzers, if scheduled.
         analyzerManager.sample(iStep_);

         // Write restart file, if scheduled.
         if (saveInterval() > 0) {
            if (iStep_ % saveInterval() == 0) {
               if (iStep_ > beginStep) {
                  simulation().save(saveFileName());
               }
            }
         }
         timer().stamp(ANALYZER);
 
         #ifdef DDMD_MODIFIERS 
         modifierManager.preIntegrate1(iStep_);
         timer().stamp(MODIFIER);
         #endif
   
         // First step of integration: Update positions, half velocity 
         integrateStep1();
         timer().stamp(INTEGRATE1);
  
         // Unset precomputed values of all energies, stresses, etc.
         simulation().modifySignal().notify();

         #ifdef DDMD_INTEGRATOR_DEBUG
         // Sanity check
         simulation().isValid();
         timer().stamp(DEBUG);
         #endif

         #ifdef DDMD_MODIFIERS 
         modifierManager.postIntegrate1(iStep_);
         timer().stamp(MODIFIER);
         #endif
   
         // Check if exchange and reneighboring is necessary
         // Note: Integrate::isExchangeNeeded uses timer.
         needExchange = isExchangeNeeded(pairPotential().skin());

         if (!atomStorage().isCartesian()) {
            UTIL_THROW("Error: atomic coordinates are not Cartesian");
         }

         // Exchange atoms if necessary
         if (needExchange) {

            #ifdef DDMD_MODIFIERS 
            modifierManager.preTransform(iStep_);
            timer().stamp(MODIFIER);
            #endif
      
            // Transform to scaled [0,1] coordinates
            atomStorage().clearSnapshot();
            atomStorage().transformCartToGen(boundary());
            timer().stamp(Integrator::TRANSFORM_F);

            #ifdef DDMD_MODIFIERS 
            modifierManager.preExchange(iStep_);
            timer().stamp(MODIFIER);
            #endif
      
            // Exchange atom ownership, reidentify ghosts
            exchanger().exchange();
            timer().stamp(Integrator::EXCHANGE);

            #ifdef DDMD_MODIFIERS 
            modifierManager.postExchange(iStep_);
            timer().stamp(MODIFIER);
            #endif
   
            // Build cell list 
            pairPotential().buildCellList();
            timer().stamp(Integrator::CELLLIST);

            // Transform from scaled [0,1] to Cartesian coordinates.
            atomStorage().transformGenToCart(boundary());
            timer().stamp(Integrator::TRANSFORM_R);

            // Build pair list
            atomStorage().makeSnapshot();
            pairPotential().buildPairList();
            timer().stamp(Integrator::PAIRLIST);

            #ifdef DDMD_MODIFIERS 
            modifierManager.postNeighbor(iStep_);
            timer().stamp(MODIFIER);
            #endif
   
         } else { // Update step (no exchange)

            #ifdef DDMD_MODIFIERS 
            modifierManager.preUpdate(iStep_);
            timer().stamp(MODIFIER);
            #endif
     
            // Update all ghost atom positions 
            exchanger().update();
            timer().stamp(UPDATE);

            #ifdef DDMD_MODIFIERS 
            modifierManager.postUpdate(iStep_);
            timer().stamp(MODIFIER);
            #endif
   
         }
         simulation().exchangeSignal().notify();
   
         #ifdef DDMD_INTEGRATOR_DEBUG
         // Sanity check
         simulation().isValid();
         timer().stamp(DEBUG);
         #endif

         if (!atomStorage().isCartesian()) {
            UTIL_THROW("Error: Atomic coordinates are not Cartesian");
         }

         #ifdef DDMD_MODIFIERS 
         modifierManager.preForce(iStep_);
         timer().stamp(MODIFIER);
         #endif
  
         // Calculate forces: 
         // Calculate new forces for all local atoms. If constant pressure
         // ensemble (not rigid), also calculate the virial stress. Both 
         // methods use the timer() internall, and both send the modifyForce 
         // signal. 
         if (simulation().boundaryEnsemble().isRigid()) {
            computeForces();
         } else {
            computeForcesAndVirial();
         }

         #ifdef DDMD_MODIFIERS 
         modifierManager.postForce(iStep_);
         timer().stamp(MODIFIER);
         #endif
   
         // 2nd step of velocity-Verlet integration. This finishes the velocity 
         // update, and normally calls simulation().velocitySignal().notify()
         integrateStep2();
         timer().stamp(INTEGRATE2);
   
         #ifdef DDMD_INTEGRATOR_DEBUG
         // Sanity check
         simulation().isValid();
         timer().stamp(DEBUG);
         #endif

         #ifdef DDMD_MODIFIERS 
         modifierManager.endOfStep(iStep_);
         timer().stamp(MODIFIER);
         #endif

      }
      exchanger().timer().stop();
      timer().stop();

      // Final analyzers and restart file, if scheduled.
      analyzerManager.sample(iStep_);
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
