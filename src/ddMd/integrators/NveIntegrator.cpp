#ifndef DDMD_NVE_INTEGRATOR_CPP
#define DDMD_NVE_INTEGRATOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "NveIntegrator.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/potentials/pair/PairPotential.h>
#include <util/space/Vector.h>
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
   * Read time step dt.
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
      pairPotential().buildCellList();
      if (!UTIL_ORTHOGONAL) {
         atomStorage().transformGenToCart(boundary());
      }
      atomStorage().makeSnapshot();
      pairPotential().buildPairList();
      simulation().computeForces();

      simulation().diagnosticManager().setup();

      // Set prefactors for acceleration
      double dtHalf = 0.5*dt_;
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = dtHalf/mass;
      }
   }

   /*
   * Integrate.
   * 
   * User must call setup() before run.
   */
   void NveIntegrator::run(int nStep)
   {
      //setup();

      Vector        dv;
      Vector        dr;
      double        prefactor; // = 0.5*dt/mass
      AtomIterator  atomIter;
      bool          needExchange;

      // Main MD loop
      nStep_ = nStep;
      timer().start();
      exchanger().timer().start();
      for (iStep_ = 0; iStep_ < nStep; ++iStep_) {

         if (Diagnostic::baseInterval > 0) {
            if (iStep_ % Diagnostic::baseInterval == 0) {
               simulation().diagnosticManager().sample(iStep_);
            }
         }
         timer().stamp(Integrator::DIAGNOSTIC);

         // 1st half of velocity Verlet.
         atomStorage().begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            prefactor = prefactors_[atomIter->typeId()];
   
            dv.multiply(atomIter->force(), prefactor);
            atomIter->velocity() += dv;
   
            dr.multiply(atomIter->velocity(), dt_);
            atomIter->position() += dr;
    
         }
         timer().stamp(Integrator::INTEGRATE1);
   
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
            timer().stamp(Integrator::UPDATE);
         }
   
         // Calculate new forces for all local atoms
         computeForces();
   
         // 2nd half of velocity Verlet
         atomStorage().begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            prefactor = prefactors_[atomIter->typeId()];
            dv.multiply(atomIter->force(), prefactor);
            atomIter->velocity() += dv;
         }
         timer().stamp(Integrator::INTEGRATE2);
   
      }
      timer().stop();
      exchanger().timer().stop();

      // Compute and reduce statistics for run.
      #ifdef UTIL_MPI
      timer().reduce(domain().communicator());
      exchanger().timer().reduce(domain().communicator());
      pairPotential().pairList().computeStatistics(domain().communicator());
      atomStorage().computeNAtomTotal(domain().communicator());
      #else
      pairPotential().pairList().computeStatistics();
      #endif

   }

}
#endif
