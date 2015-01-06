/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
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
     vcoeff_(),
     fcoeff_(),
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
      vcoeff_.allocate(nAtomType);
      fcoeff_.allocate(nAtomType);
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
      ar & vcoeff_;
      ar & fcoeff_;
   }

   /*
   * Save the internal state to an archive.
   */
   void NvtLangevinIntegrator::save(Serializable::OArchive& ar)
   {
      ar & dt_;
      ar & gamma_;
      ar & prefactors_;
      ar & vcoeff_;
      ar & fcoeff_;
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
      double temp = energyEnsemble.temperature();
      double vc = (exp(-dt_*gamma_) - 1.0)/dt_;
      double fc = 12.0*temp*(1.0 - exp(-2.0*dt_*gamma_))/(dt_*dt_);

      // Loop over atom types
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = 0.5*dt_/mass;
         vcoeff_[i] = mass*vc;
         fcoeff_[i] = sqrt(mass*fc);
      }

   }

   /*
   * Verlet MD NVE integrator step
   *
   * This method implements the algorithm:
   *
   *        vm(n)  = v(n) + 0.5*a(n)*dt
   *        x(n+1) = x(n) + vm(n)*dt
   *
   *        calculate determinstic force a(n+1)
   *        add Langevin force
   *      
   *        v(n+1) = vm(n) + 0.5*a(n+1)*dt
   *
   * where x is position, v is velocity, and a is acceleration.
   */
   void NvtLangevinIntegrator::step()
   {
      Vector dr;
      Vector dv;
      Vector df;
      System::MoleculeIterator molIter;
      double fc;
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

               // Velocity update (half step)
               dv.multiply(atomIter->force(), prefactors_[typeId]);
               atomIter->velocity() += dv;

               // Position update (full step)
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

               // Velocity update (half step)
               dv.multiply(atomPtr->force(), prefactors_[typeId]);
               atomPtr->velocity() += dv;

               // Position update (full step)
               dr.multiply(atomPtr->velocity(), dt_);
               atomPtr->position() += dr;
            }
         }
      }
      #endif

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

               // Add Langevin force to atomic force
               df.multiply(atomIter->velocity(), vcoeff_[typeId]);
               fc = fcoeff_[typeId];
               for (j=0; j < Dimension; ++j) {
                  df[j] += (random.uniform() - 0.5)*fc;
               }
               atomIter->force() += df;

               // Velocity update
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

               // Add to Langevin force to atomic force
               df.multiply(atomPtr->velocity(), vcoeff_[typeId]);
               fc = fcoeff_[typeId];
               for (j=0; j < Dimension; ++j) {
                  df[j] += (random.uniform() - 0.5)*fc;
               }
               atomPtr->force() += df;

               // Velocity update
               dv.multiply(atomPtr->force(), prefactors_[typeId]);
               atomPtr->velocity() += dv;
            }
         }
      }
      #endif

      #ifndef INTER_NOPAIR
      if (!system().pairPotential().isPairListCurrent()) {
         system().pairPotential().buildPairList();
      }
      #endif

   }

}
#undef USE_ITERATOR
