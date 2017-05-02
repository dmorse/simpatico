/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NvtLangevinIntegrator.h"
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <util/boundary/Boundary.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/archives/Serializable_includes.h>

//#define USE_ITERATOR 1
namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   NvtLangevinIntegrator::NvtLangevinIntegrator(MdSystem& system)
   : MdIntegrator(system),
     prefactors_(),
     cv_(),
     cr_(),
     gamma_()
   {  setClassName("NvtLangevinIntegrator"); }

   /*
   * Destructor.
   */
   NvtLangevinIntegrator::~NvtLangevinIntegrator()
   {}

   /*
   * Read parameter and configuration files, initialize system.
   */
   void NvtLangevinIntegrator::readParameters(std::istream &in)
   {
      read<double>(in, "dt", dt_);
      read<double>(in, "gamma", gamma_);
      int nAtomType = simulation().nAtomType();
      prefactors_.allocate(nAtomType);
      cv_.allocate(nAtomType);
      cr_.allocate(nAtomType);
   }

   /*
   * Load the internal state to an archive.
   */
   void NvtLangevinIntegrator::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<double>(ar, "dt", dt_);
      int nAtomType = simulation().nAtomType();
      prefactors_.allocate(nAtomType);
      ar & prefactors_;
      ar & cv_;
      ar & cr_;
   }

   /*
   * Save the internal state to an archive.
   */
   void NvtLangevinIntegrator::save(Serializable::OArchive& ar)
   {
      ar & dt_;
      ar & gamma_;
      ar & prefactors_;
      ar & cv_;
      ar & cr_;
   }

   /*
   * Initialize constants.
   */
   void NvtLangevinIntegrator::setup()
   {
      const EnergyEnsemble& energyEnsemble = system().energyEnsemble();
      if (!energyEnsemble.isIsothermal()) {
         UTIL_THROW("Energy ensemble is not isothermal");
      }
      double cv = (exp(-dt_*gamma_) - 1.0)/dt_;
      double temp = energyEnsemble.temperature();
      double d = 2.0/(1.0 + exp(-dt_*gamma_));
      double cr = 12.0*temp*d*(1.0 - exp(-2.0*dt_*gamma_))/(dt_*dt_);

      // Loop over atom types
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = 0.5*dt_/mass;
         cv_[i] = mass*cv;
         cr_[i] = sqrt(mass*cr);
      }
      system().positionSignal().notify();
      system().velocitySignal().notify();

   }

   /*
   * Verlet MD NVE integrator step
   *
   * This method implements the algorithm:
   *
   *        vm(n)  = v(n) + 0.5*a(n)*dt
   *        x(n+1) = x(n) + vm(n)*dt
   *
   *        calculate determinstic force f(n+1)
   *        add Langevin force, calculated using vm(n)
   *      
   *        v(n+1) = vm(n) + 0.5*f(n+1)*dt/m
   *
   * where x is position, v is velocity, and f is force.
   */
   void NvtLangevinIntegrator::step()
   {
      Vector dr;
      Vector dv;
      Vector df;
      double cr;
      System::MoleculeIterator molIter;
      int iSpecies, nSpecies, typeId;

      nSpecies = simulation().nSpecies();

      // 1st half velocity Verlet, loop over atoms
      #if USE_ITERATOR
      Molecule::AtomIterator atomIter;
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            molIter->begin(atomIter); 
            for ( ; atomIter.notEnd(); ++atomIter) {
               typeId = atomIter->typeId();

               // Update velocity one half step, prefactor = dt/(2m)
               dv.multiply(atomIter->force(), prefactors_[typeId]);
               atomIter->velocity() += dv;

               // Update position (full step)
               dr.multiply(atomIter->velocity(), dt_);
               atomIter->position() += dr;
            }
         }
      }
      #else
      Atom* atomPtr;
      int ia;
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            for (ia=0; ia < molIter->nAtom(); ++ia) {
               atomPtr = &molIter->atom(ia);
               typeId = atomPtr->typeId();

               // Update velocity one half step, prefactor = dt/(2m)
               dv.multiply(atomPtr->force(), prefactors_[typeId]);
               atomPtr->velocity() += dv;

               // Position update (full step)
               dr.multiply(atomPtr->velocity(), dt_);
               atomPtr->position() += dr;
            }
         }
      }
      #endif
      system().positionSignal().notify();
      system().velocitySignal().notify();

      #ifndef SIMP_NOPAIR
      if (!system().pairPotential().isPairListCurrent()) {
         system().pairPotential().buildPairList();
      }
      #endif

      // Calculate conservative force
      system().calculateForces();

      // 2nd half velocity Verlet, loop over atoms
      Random& random = simulation().random();
      int j;
      #if USE_ITERATOR
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            molIter->begin(atomIter); 
            for ( ; atomIter.notEnd(); ++atomIter) {
               typeId = atomIter->typeId();

               // Add Langevin drag and random force to atomic force
               df.multiply(atomIter->velocity(), cv_[typeId]);
               cr = cr_[typeId];
               for (j=0; j < Dimension; ++j) {
                  df[j] += (random.uniform() - 0.5)*cr;
               }
               atomIter->force() += df;

               // Velocity update by half step
               dv.multiply(atomIter->force(), prefactors_[typeId]);
               atomIter->velocity() += dv;
            }
         }
      }
      #else
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            for (ia=0; ia < molIter->nAtom(); ++ia) {
               atomPtr = &molIter->atom(ia);
               typeId = atomPtr->typeId();

               // Add Langevin drag and random force to atomic force
               df.multiply(atomPtr->velocity(), cv_[typeId]);
               cr = cr_[typeId];
               for (j=0; j < Dimension; ++j) {
                  df[j] += (random.uniform() - 0.5)*cr;
               }
               atomPtr->force() += df;

               // Velocity update by half step
               dv.multiply(atomPtr->force(), prefactors_[typeId]);
               atomPtr->velocity() += dv;
            }
         }
      }
      #endif
      system().velocitySignal().notify();

   }

}
#undef USE_ITERATOR
