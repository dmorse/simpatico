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
#include <ddMd/interaction/Interaction.h>
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
      Interaction* interactionPtr  = &systemPtr_->interaction();

      atomStoragePtr->clearSnapshot();
      //atomStoragePtr->clearGhosts();
      exchangerPtr->exchangeAtoms();
      exchangerPtr->exchangeGhosts();
      atomStoragePtr->makeSnapshot();
      interactionPtr->findNeighbors();
      interactionPtr->calculateForces();
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

      AtomStorage* atomStoragePtr  = &systemPtr_->atomStorage();
      Exchanger*   exchangerPtr    = &systemPtr_->exchanger();
      Interaction* interactionPtr  = &systemPtr_->interaction();

      // 1st half of velocity Verlet.
      atomStoragePtr->begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {

         dv.multiply(atomIter->force(), dtHalf);
         atomIter->velocity() += dv;

         dr.multiply(atomIter->velocity(), dt_);
         atomIter->position() += dr;
 
      }

      bool needExchange = systemPtr_->needExchange();


      // Exchange atoms if necessary
      if (needExchange) {

         atomStoragePtr->clearSnapshot();
         //atomStoragePtr->clearGhosts();
         exchangerPtr->exchangeAtoms();
         exchangerPtr->exchangeGhosts();
         atomStoragePtr->makeSnapshot();
         interactionPtr->findNeighbors();

      } else {

         exchangerPtr->updateGhosts();

      }

      // Calculate new forces for all local atoms
      interactionPtr->calculateForces();

      // 2nd half of velocity Verlet
      atomStoragePtr->begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {

         dv.multiply(atomIter->force(), dtHalf);
         atomIter->velocity() += dv;

      }

   }

}
#endif
