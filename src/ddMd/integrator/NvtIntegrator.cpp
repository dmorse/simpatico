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
#include <ddMd/system/System.h>
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
   NvtIntegrator::NvtIntegrator(System& system)
    : Integrator(system),
      prefactors_(),
      T_target_(1.0),
      T_kinetic_(1.0),
      xi_(0.0),
      xiDot_(0.0),
      tauT_(1.0),
      nuT_(1.0)
   {
      // Note: Within the constructor, the method parameter "system" hides 
      // the system() method name.

      // Precondition
      if (!system.energyEnsemble().isIsothermal() ) {
         UTIL_THROW("System energy ensemble is not isothermal");
      }
   }

   /*
   * Destructor.
   */
   NvtIntegrator::~NvtIntegrator()
   {}

   /* 
   * Read parameter and configuration files, initialize system.
   */
   void NvtIntegrator::readParam(std::istream &in) 
   {
      //readBegin(in, "NvtIntegrator");
      read<double>(in, "dt",   dt_);
      read<double>(in, "tauT", tauT_);
      //readEnd(in);

      nuT_ = 1.0/tauT_;
      int nAtomType = system().nAtomType();
      if (!prefactors_.isAllocated()) {
         prefactors_.allocate(nAtomType);
      }

   }

   void NvtIntegrator::setup()
   {
      AtomStorage* atomStoragePtr = &system().atomStorage();
      Domain*      domainPtr = &system().domain();
      Exchanger*   exchangerPtr = &system().exchanger();
      PairPotential* pairPotentialPtr = &system().pairPotential();


      double dtHalf = 0.5*dt_;
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = system().atomType(i).mass();
         prefactors_[i] = dtHalf/mass;
      }

      atomStoragePtr->clearSnapshot();
      exchangerPtr->exchange();
      atomStoragePtr->makeSnapshot();
      pairPotentialPtr->findNeighbors();
      system().computeForces();

      // Initialize nAtom_, xiDot_, xi_
      system().computeKineticEnergy();
      atomStoragePtr->computeNAtomTotal(domainPtr->communicator());
      if (domainPtr->isMaster()) {
         T_target_ = system().energyEnsemble().temperature();
         nAtom_  = atomStoragePtr->nAtomTotal();
         T_kinetic_ = system().kineticEnergy()*2.0/double(3*nAtom_);
         xiDot_ = (T_kinetic_/T_target_ -1.0)*nuT_*nuT_;
      }
      bcast(domainPtr->communicator(), xiDot_, 0);
      xi_ = 0.0;
   }

   /*
   * Nose-Hoover integrator step.
   *
   * This implements a reversible Velocity-Verlet MD NVT integrator step.
   *
   * Reference: Winkler, Kraus, and Reineker, J. Chem. Phys. 102, 9018 (1995).
   */
   void NvtIntegrator::step() 
   {
      Vector  dv;
      Vector  dr;
      double  dtHalf = 0.5*dt_;
      double  prefactor;
      double  factor;
      AtomIterator atomIter;

      Domain*  domainPtr = &system().domain();
      AtomStorage* atomStoragePtr = &system().atomStorage();
      Exchanger* exchangerPtr = &system().exchanger();
      PairPotential* pairPotentialPtr = &system().pairPotential();

      T_target_ = system().energyEnsemble().temperature();
      factor = exp(-dtHalf*(xi_ + xiDot_*dtHalf));

      // 1st half of velocity Verlet.
      atomStoragePtr->begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         atomIter->velocity() *= factor;
         prefactor = prefactors_[atomIter->typeId()];
         dv.multiply(atomIter->force(), prefactor);
         atomIter->velocity() += dv;
         dr.multiply(atomIter->velocity(), dt_);
         atomIter->position() += dr;
      }

      // Exchange atoms if necessary
      if (system().needExchange()) {
         atomStoragePtr->clearSnapshot();
         exchangerPtr->exchange();
         atomStoragePtr->makeSnapshot();
         pairPotentialPtr->findNeighbors();
      } else {
         exchangerPtr->update();
      }

      // Calculate new forces for all local atoms
      system().computeForces();

      // 2nd half of velocity Verlet
      atomStoragePtr->begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         prefactor = prefactors_[atomIter->typeId()];
         dv.multiply(atomIter->force(), prefactor);
         atomIter->velocity() += dv;
         atomIter->velocity() *=factor;
      }

      // Update xiDot_ and xi_
      system().computeKineticEnergy();
      if (domainPtr->isMaster()) {
         xi_ += xiDot_*dtHalf;
         T_kinetic_ = system().kineticEnergy()*2.0/double(3*nAtom_);
         xiDot_ = (T_kinetic_/T_target_  - 1.0)*nuT_*nuT_;
         xi_ += xiDot_*dtHalf;
      }
      bcast(domainPtr->communicator(), xiDot_, 0);
      bcast(domainPtr->communicator(), xi_, 0);
   }

}
#endif
