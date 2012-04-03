#ifndef DDMD_NVE_INTEGRATOR_CPP
#define DDMD_NVE_INTEGRATOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "NveIntegrator.h"
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
   NveIntegrator::NveIntegrator(Simulation& simulation)
    : Integrator(simulation)
   {}

   /*
   * Destructor.
   */
   NveIntegrator::~NveIntegrator()
   {}

   /*
   * Read time step.
   */
   void NveIntegrator::readParam(std::istream& in)
   {
      //readBegin(in, "NveIntegrator");
      read<double>(in, "dt", dt_);
      //readEnd(in);

      int nAtomType = simulation().nAtomType();
      if (!prefactors_.isAllocated()) {
         prefactors_.allocate(nAtomType);
      }

   }

   void NveIntegrator::setup()
   {
      atomStorage().clearSnapshot();
      exchanger().exchange();
      atomStorage().makeSnapshot();
      pairPotential().findNeighbors();
      simulation().computeForces();

      double dtHalf = 0.5*dt_;
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = dtHalf/mass;
      }
   }

   void NveIntegrator::step()
   {
      // Preconditions
      //if (!storage.isInitialized()) {
      //   UTIL_THROW("AtomStorage must be initialized");
      //}

      Vector        dv;
      Vector        dr;
      double        prefactor; // = 0.5*dt/mass
      AtomIterator  atomIter;

      // 1st half of velocity Verlet.
      atomStorage().begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         prefactor = prefactors_[atomIter->typeId()];

         dv.multiply(atomIter->force(), prefactor);
         atomIter->velocity() += dv;

         dr.multiply(atomIter->velocity(), dt_);
         atomIter->position() += dr;
      }

      // Exchange atoms if necessary
      if (simulation().needExchange()) {
         atomStorage().clearSnapshot();
         exchanger().exchange();
         atomStorage().makeSnapshot();
         pairPotential().findNeighbors();
      } else {
         exchanger().update();
      }

      // Calculate new forces for all local atoms
      simulation().computeForces();

      // 2nd half of velocity Verlet
      atomStorage().begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         prefactor = prefactors_[atomIter->typeId()];
         dv.multiply(atomIter->force(), prefactor);
         atomIter->velocity() += dv;
      }

   }

   /*
   * Integrate.
   */
   void NveIntegrator::run(int nStep)
   {
      if(domain().isMaster()) {
         Log::file() << std::endl;
      }
      nStep_ = nStep;

      // Set prefactor_[i] = 0.5*dt/mass for each atom type i.
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = 0.5*dt_/mass;
      }

      atomStorage().clearSnapshot();
      exchanger().exchange();
      atomStorage().makeSnapshot();
      pairPotential().findNeighbors();
      simulation().computeForces();

      simulation().diagnosticManager().setup();

      Vector        dv;
      Vector        dr;
      double        prefactor; // = 0.5*dt/mass
      AtomIterator  atomIter;

      // Main MD loop
      timer().start();
      for (iStep_ = 0; iStep_ < nStep; ++iStep_) {

         if (Diagnostic::baseInterval > 0) {
            if (iStep_ % Diagnostic::baseInterval == 0) {
               simulation().diagnosticManager().sample(iStep_);
            }
         }

         // 1st half of velocity Verlet.
         atomStorage().begin(atomIter);
         for ( ; !atomIter.atEnd(); ++atomIter) {
            prefactor = prefactors_[atomIter->typeId()];
   
            dv.multiply(atomIter->force(), prefactor);
            atomIter->velocity() += dv;
   
            dr.multiply(atomIter->velocity(), dt_);
            atomIter->position() += dr;
    
         }
   
         // Exchange atoms if necessary
         if (simulation().needExchange()) {
            atomStorage().clearSnapshot();
            exchanger().exchange();
            atomStorage().makeSnapshot();
            pairPotential().findNeighbors();
         } else {
            exchanger().update();
         }
   
         // Calculate new forces for all local atoms
         simulation().computeForces();
   
         // 2nd half of velocity Verlet
         atomStorage().begin(atomIter);
         for ( ; !atomIter.atEnd(); ++atomIter) {
            prefactor = prefactors_[atomIter->typeId()];
            dv.multiply(atomIter->force(), prefactor);
            atomIter->velocity() += dv;
         }
   
      }
      timer().stop();
   }

}
#endif
