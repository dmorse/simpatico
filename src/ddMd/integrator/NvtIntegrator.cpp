#ifndef DDMD_NVT_INTEGRATOR_CPP
#define DDMD_NVT_INTEGRATOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "NvtIntegrator.h"
#include <util/ensembles/EnergyEnsemble.h>
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/potentials/pair/PairPotential.h>
#include <util/space/Vector.h>
#include <util/util/Timer.h>
#include <util/global.h>

#include <iostream>

namespace DdMd
{
   using namespace Util;

   /* 
   * Constructor.
   */
   NvtIntegrator::NvtIntegrator(Simulation& simulation)
    : Integrator(simulation),
      prefactors_(),
      T_target_(1.0),
      T_kinetic_(1.0),
      xi_(0.0),
      xiDot_(0.0),
      tauT_(1.0),
      nuT_(1.0)
   {
      // Note: Within the constructor, the method parameter "simulation" hides 
      // the simulation() method name.

      // Precondition
      if (!simulation.energyEnsemble().isIsothermal() ) {
         UTIL_THROW("Simulation energy ensemble is not isothermal");
      }
   }

   /*
   * Destructor.
   */
   NvtIntegrator::~NvtIntegrator()
   {}

   /* 
   * Read parameter and configuration files, initialize simulation.
   */
   void NvtIntegrator::readParam(std::istream &in) 
   {
      //readBegin(in, "NvtIntegrator");
      read<double>(in, "dt",   dt_);
      read<double>(in, "tauT", tauT_);
      //readEnd(in);

      nuT_ = 1.0/tauT_;
      int nAtomType = simulation().nAtomType();
      if (!prefactors_.isAllocated()) {
         prefactors_.allocate(nAtomType);
      }

   }

   void NvtIntegrator::setup()
   {
      double dtHalf = 0.5*dt_;
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = dtHalf/mass;
      }

      atomStorage().clearSnapshot();
      exchanger().exchange();
      atomStorage().makeSnapshot();
      pairPotential().findNeighbors();
      simulation().computeForces();

      // Initialize nAtom_, xiDot_, xi_
      simulation().computeKineticEnergy();
      #ifdef UTIL_MPI
      atomStorage().computeNAtomTotal(domain().communicator());
      #endif
      if (domain().isMaster()) {
         T_target_ = simulation().energyEnsemble().temperature();
         nAtom_  = atomStorage().nAtomTotal();
         T_kinetic_ = simulation().kineticEnergy()*2.0/double(3*nAtom_);
         xiDot_ = (T_kinetic_/T_target_ -1.0)*nuT_*nuT_;
      }
      #ifdef UTIL_MPI
      bcast(domain().communicator(), xiDot_, 0);
      #endif
      xi_ = 0.0;
   }

   /*
   * Integrate Nose-Hoover.
   *
   * This implements a reversible Velocity-Verlet MD NVT integrator step.
   *
   * Reference: Winkler, Kraus, and Reineker, J. Chem. Phys. 102, 9018 (1995).
   */
   void NvtIntegrator::run(int nStep)
   {
      nStep_ = nStep;
      if (domain().isMaster()) {
         Log::file() << std::endl;
      }

      setup();
      simulation().diagnosticManager().setup();

      Vector dv;
      Vector dr;
      double prefactor; // = 0.5*dt/mass
      double dtHalf = 0.5*dt_;
      double factor;
      AtomIterator atomIter;
      bool needExchange;


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
   
         T_target_ = simulation().energyEnsemble().temperature();
         factor = exp(-dtHalf*(xi_ + xiDot_*dtHalf));
   
         // 1st half of velocity Verlet.
         atomStorage().begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            atomIter->velocity() *= factor;
            prefactor = prefactors_[atomIter->typeId()];
            dv.multiply(atomIter->force(), prefactor);
            atomIter->velocity() += dv;
            dr.multiply(atomIter->velocity(), dt_);
            atomIter->position() += dr;
         }
         timer().stamp(INTEGRATE1);
   
   
         // Check if exchange and reneighboring is necessary
         needExchange = simulation().needExchange();
         timer().stamp(Integrator::CHECK);

         // Exchange atoms if necessary
         if (needExchange) {
            atomStorage().clearSnapshot();
            exchanger().exchange();
            timer().stamp(EXCHANGE);
            atomStorage().makeSnapshot();
            pairPotential().findNeighbors();
            timer().stamp(NEIGHBOR);
         } else {
            exchanger().update();
            timer().stamp(UPDATE);
         }
   
         // Calculate new forces for all local atoms
         simulation().computeForces();
         timer().stamp(FORCE);
   
         // 2nd half of velocity Verlet
         atomStorage().begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            prefactor = prefactors_[atomIter->typeId()];
            dv.multiply(atomIter->force(), prefactor);
            atomIter->velocity() += dv;
            atomIter->velocity() *=factor;
         }
   
         // Update xiDot_ and xi_
         simulation().computeKineticEnergy();
         if (domain().isMaster()) {
            xi_ += xiDot_*dtHalf;
            T_kinetic_ = simulation().kineticEnergy()*2.0/double(3*nAtom_);
            xiDot_ = (T_kinetic_/T_target_  - 1.0)*nuT_*nuT_;
            xi_ += xiDot_*dtHalf;
         }
         #ifdef UTIL_MPI
         bcast(domain().communicator(), xiDot_, 0);
         bcast(domain().communicator(), xi_, 0);
         #endif
         timer().stamp(INTEGRATE2);
   
      }
      exchanger().timer().stop();
      timer().stop();

   }

}
#endif
