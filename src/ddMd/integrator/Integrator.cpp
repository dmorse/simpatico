#ifndef INTEGRATOR_CPP
#define INTEGRATOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Integrator.h"
#include <ddMd/system/System.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/potentials/PairPotential.h>
#include <util/space/Vector.h>
#include <util/global.h>

#include <iostream>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   Integrator::Integrator(System& system)
    : systemPtr_(&system)
   {}

   /*
   * Destructor.
   */
   Integrator::~Integrator()
   {}

   /*
   * Read time step.
   */
   void Integrator::readParam(std::istream& in)
   {
      readBegin(in, "Integrator");
      read<double>(in, "dt", dt_);
      readEnd(in);
   }

   void Integrator::initialize()
   {
      AtomStorage* atomStoragePtr  = &systemPtr_->atomStorage();
      Exchanger*   exchangerPtr    = &systemPtr_->exchanger();
      PairPotential* pairPotentialPtr  = &systemPtr_->pairPotential();

      atomStoragePtr->clearSnapshot();
      exchangerPtr->exchange();
      atomStoragePtr->makeSnapshot();
      pairPotentialPtr->findNeighbors();
      systemPtr_->computeForces();
   }

   void Integrator::step()
   {
      // Preconditions
      //if (!storage.isInitialized()) {
      //   UTIL_THROW("AtomStorage must be initialized");
      //}

      Vector        dv;
      Vector        dr;
      double        dtHalf = 0.5*dt_;
      AtomIterator  atomIter;

      AtomStorage* atomStoragePtr = &systemPtr_->atomStorage();
      Exchanger* exchangerPtr = &systemPtr_->exchanger();
      PairPotential* pairPotentialPtr = &systemPtr_->pairPotential();

      // 1st half of velocity Verlet.
      atomStoragePtr->begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {

         dv.multiply(atomIter->force(), dtHalf);
         atomIter->velocity() += dv;

         dr.multiply(atomIter->velocity(), dt_);
         atomIter->position() += dr;
 
      }

      // Exchange atoms if necessary
      if (systemPtr_->needExchange()) {
         atomStoragePtr->clearSnapshot();
         exchangerPtr->exchange();
         atomStoragePtr->makeSnapshot();
         pairPotentialPtr->findNeighbors();
      } else {
         exchangerPtr->update();
      }

      // Calculate new forces for all local atoms
      systemPtr_->computeForces();

      // 2nd half of velocity Verlet
      atomStoragePtr->begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         dv.multiply(atomIter->force(), dtHalf);
         atomIter->velocity() += dv;
      }

   }

}
#endif
