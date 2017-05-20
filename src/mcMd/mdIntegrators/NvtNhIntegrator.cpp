/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NvtNhIntegrator.h"
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/archives/Serializable_includes.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   NvtNhIntegrator::NvtNhIntegrator(MdSystem& system)
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

      setClassName("NvtNhIntegrator"); 
   }

   /* 
   * Destructor.   
   */
   NvtNhIntegrator::~NvtNhIntegrator() 
   {}

   /* 
   * Read parameter and configuration files, initialize system.
   */
   void NvtNhIntegrator::readParameters(std::istream &in) 
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

   /*
   * Load the internal state to an archive.
   */
   void NvtNhIntegrator::loadParameters(Serializable::IArchive& ar)
   {  
      loadParameter<double>(ar, "dt",   dt_);
      loadParameter<double>(ar, "tauT", tauT_);
      ar & nuT_;
      ar & T_target_;
      ar & T_kinetic_;
      ar & xiDot_;
      ar & xi_;
      ar & prefactors_;
   }

   /*
   * Save the internal state to an archive.
   */
   void NvtNhIntegrator::save(Serializable::OArchive& ar)
   {
      ar & dt_;
      ar & tauT_;
      ar & nuT_;
      ar & T_target_;
      ar & T_kinetic_;
      ar & xiDot_;
      ar & xi_;
      ar & prefactors_;
   }

   void NvtNhIntegrator::setup() 
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
      system().positionSignal().notify();
      system().velocitySignal().notify();
   }

   /*
   * Nose-Hoover integrator step.
   *
   * This implements a reversible Velocity-Verlet MD NVT integrator step.
   *
   * Reference: Winkler, Kraus, and Reineker, J. Chem. Phys. 102, 9018 (1995).
   */
   void NvtNhIntegrator::step() 
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
      system().positionSignal().notify();
      system().velocitySignal().notify();

      // First half of update of xi_
      xi_ += xiDot_*dtHalf;

      #ifndef SIMP_NOPAIR
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
      system().velocitySignal().notify();

      // Update xiDot and complete update of xi_
      T_kinetic_ = system().kineticEnergy()*2.0/double(3*nAtom);
      xiDot_ = (T_kinetic_/T_target_ -1.0)*nuT_*nuT_;
      xi_ += xiDot_*dtHalf;

   }

}
