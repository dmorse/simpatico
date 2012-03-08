#ifndef NVE_INTEGRATOR_CPP
#define NVE_INTEGRATOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "NvtIntegrator.h"
#include <ddMd/ensembles/EnergyEnsemble.h>
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
    : Integrator(system)
   {}

   /* 
   * Constructor.
   */
   NvtIntegrator::NvtIntegrator(MdSystem& system)
    : MdIntegrator(system),
      T_target_(1.0),
      T_kinetic_(1.0),
      xi_(0.0),
      xiDot_(0.0),
      tauT_(1.0),
      nuT_(1.0),
      energyEnsemblePtr_(0)
   {
      // Note: Within the constructor, the method parameter "system" hides 
      // the system() method name.

      // Precondition
      if (!system.energyEnsemble().isIsothermal() ) {
         UTIL_THROW("System energy ensemble is not isothermal");
      }

      energyEnsemblePtr_ = &(system.energyEnsemble());
      T_target_          = energyEnsemblePtr_->temperature();
      T_kinetic_         = T_target_;

   }

   /*
   * Destructor.
   */
   NvtIntegrator::~NvtIntegrator()
   {}

   /*
   * Read time step.
   */
   void NvtIntegrator::readParam(std::istream& in)
   {
      readBegin(in, "NvtIntegrator");
      read<double>(in, "dt", dt_);
      readEnd(in);

      int nAtomType = system().nAtomType();
      if (!prefactors_.isAllocated()) {
         prefactors_.allocate(nAtomType);
      }

   }
   /* 
   * Read parameter and configuration files, initialize system.
   */
   void NvtIntegrator::readParam(std::istream &in) 
   {
      read<double>(in, "dt",   dt_);
      read<double>(in, "tauT", tauT_);
      nuT_ = 1.0/tauT_;
      T_target_  = energyEnsemblePtr_->temperature();
      T_kinetic_ = T_target_;
      xiDot_ = 0.0;
      xi_    = 0.0;

      int nAtomType = simulation().nAtomType();
      if (!prefactors_.isAllocated()) {
         prefactors_.allocate(nAtomType);
      }

   }

   void NvtIntegrator::setup()
   {
      AtomStorage* atomStoragePtr = &system().atomStorage();
      Exchanger*   exchangerPtr = &system().exchanger();
      PairPotential* pairPotentialPtr = &system().pairPotential();

      atomStoragePtr->clearSnapshot();
      exchangerPtr->exchange();
      atomStoragePtr->makeSnapshot();
      pairPotentialPtr->findNeighbors();
      system().computeForces();

      double dtHalf = 0.5*dt_;
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = system().atomType(i).mass();
         prefactors_[i] = dtHalf/mass;
      }
   }

   void NvtIntegrator::setup() 
   {
      double mass, dtHalf;
      int nAtomType = simulation().nAtomType();
      int nAtom  = system().nAtom();

      T_kinetic_ = system().kineticEnergy()*2.0/double(3*nAtom);
      T_target_ = energyEnsemblePtr_->temperature();
      xiDot_ = (T_kinetic_/T_target_ -1.0)*nuT_*nuT_;
      xi_ = 0.0;

      dtHalf = 0.5*dt_;
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = dtHalf/mass;
      }
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
      System::MoleculeIterator molIter;
      double  dtHalf = 0.5*dt_;
      double  prefactor;
      double  factor;
      Molecule::AtomIterator atomIter;
      int  iSpecies, nSpecies;
      int  nAtom;

      T_target_ = energyEnsemblePtr_->temperature();
      nSpecies  = simulation().nSpecies();
      nAtom     = system().nAtom();

      factor = exp(-dtHalf*(xi_ + xiDot_*dtHalf));

      // 1st half velocity Verlet, loop over atoms 
      for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {

            molIter->begin(atomIter); 
            for ( ; atomIter.notEnd(); ++atomIter) {

               atomIter->velocity() *= factor;

               prefactor = prefactors_[atomIter->typeId()];
               dv.multiply(atomIter->force(), prefactor);
               //dv.multiply(atomIter->force(), dtHalf);

               atomIter->velocity() += dv;
               dr.multiply(atomIter->velocity(), dt_);

               atomIter->position() += dr;

            }

         }
      }

      // First half of update of xi_
      xi_ += xiDot_*dtHalf;

      #ifndef MCMD_NOPAIR
      // Rebuild the pair list if necessary
      if (!system().pairPotential().isPairListCurrent()) {
         system().pairPotential().buildPairList();
      }
      #endif

      system().calculateForces();

      // 2nd half velocity Verlet, loop over atoms
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               prefactor = prefactors_[atomIter->typeId()];
               dv.multiply(atomIter->force(), prefactor);
               atomIter->velocity() += dv;
               atomIter->velocity() *=factor;
            }
         }
      }

      // Update xiDot and complete update of xi_
      T_kinetic_ = system().kineticEnergy()*2.0/double(3*nAtom);
      xiDot_     = (T_kinetic_/T_target_ -1.0)*nuT_*nuT_;
      xi_       += xiDot_*dtHalf;

   }

   void NvtIntegrator::step()
   {
      // Preconditions
      //if (!storage.isInitialized()) {
      //   UTIL_THROW("AtomStorage must be initialized");
      //}

      /// DdMd::Nve step

      Vector        dv;
      Vector        dr;
      double        prefactor; // = 0.5*dt/mass
      AtomIterator  atomIter;

      AtomStorage* atomStoragePtr = &system().atomStorage();
      Exchanger* exchangerPtr = &system().exchanger();
      PairPotential* pairPotentialPtr = &system().pairPotential();

      // 1st half of velocity Verlet.
      atomStoragePtr->begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
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
      }

   }

}
#endif
