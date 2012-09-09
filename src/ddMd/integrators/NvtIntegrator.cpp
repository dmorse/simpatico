#ifndef DDMD_NVT_INTEGRATOR_CPP
#define DDMD_NVT_INTEGRATOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
    : TwoStepIntegrator(simulation),
      prefactors_(),
      T_target_(1.0),
      T_kinetic_(1.0),
      xi_(0.0),
      xiDot_(0.0),
      tauT_(1.0),
      nuT_(1.0)
   {
      setClassName("NvtIntegrator");

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
   void NvtIntegrator::readParameters(std::istream &in) 
   {
      read<double>(in, "dt",   dt_);
      read<double>(in, "tauT", tauT_);

      nuT_ = 1.0/tauT_;
      int nAtomType = simulation().nAtomType();
      if (!prefactors_.isAllocated()) {
         prefactors_.allocate(nAtomType);
      }

   }

   void NvtIntegrator::setup()
   {
      #if 0
      atomStorage().clearSnapshot();
      exchanger().exchange();
      pairPotential().buildCellList();
      if (!UTIL_ORTHOGONAL) {
         atomStorage().transformGenToCart(boundary());
      }
      pairPotential().buildPairList();
      atomStorage().makeSnapshot();
      simulation().computeForces();
      #endif

      // Exchange atoms, build pair list, compute forces.
      setupAtoms();

      simulation().diagnosticManager().setup();

      // Calculate prefactors for acceleration
      double dtHalf = 0.5*dt_;
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = dtHalf/mass;
      }

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

   void NvtIntegrator::integrateStep1()
   {
      Vector dv;
      Vector dr;
      double prefactor; // = 0.5*dt/mass
      double dtHalf = 0.5*dt_;
      double factor;
      AtomIterator atomIter;

   
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
   }

   void NvtIntegrator::integrateStep2()
   {
      Vector dv;
      double prefactor; // = 0.5*dt/mass
      double dtHalf = 0.5*dt_;
      double factor;
      AtomIterator atomIter;

      T_target_ = simulation().energyEnsemble().temperature();
      factor = exp(-dtHalf*(xi_ + xiDot_*dtHalf));

      // 2nd half of velocity Verlet
      atomStorage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         prefactor = prefactors_[atomIter->typeId()];
         dv.multiply(atomIter->force(), prefactor);
         atomIter->velocity() += dv;
         atomIter->velocity() *=factor;
      }

      // Notify observers of change in velocity
      simulation().velocitySignal().notify();

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

   }

}
#endif
