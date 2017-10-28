/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NvtLangevinIntegrator.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/potentials/pair/PairPotential.h>
#include <simp/ensembles/EnergyEnsemble.h>
#include <util/space/Vector.h>
#include <util/global.h>

#include <iostream>

namespace DdMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   NvtLangevinIntegrator::NvtLangevinIntegrator(Simulation& simulation)
    : TwoStepIntegrator(simulation),
     dt_(0.0),
     gamma_(0.0),
     prefactors_(),
     cv_(),
     cr_()
   {  setClassName("NvtLangevinIntegrator"); }

   /*
   * Destructor.
   */
   NvtLangevinIntegrator::~NvtLangevinIntegrator()
   {}

   /*
   * Read time step dt.
   */
   void NvtLangevinIntegrator::readParameters(std::istream& in)
   {
      read<double>(in, "dt", dt_);
      read<double>(in, "gamma", gamma_);
      Integrator::readParameters(in);

      int nAtomType = simulation().nAtomType();
      if (!prefactors_.isAllocated()) {
         prefactors_.allocate(nAtomType);
         cv_.allocate(nAtomType);
         cr_.allocate(nAtomType);
      }
   }

   /**
   * Load internal state from an archive.
   */
   void NvtLangevinIntegrator::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<double>(ar, "dt", dt_);
      Integrator::loadParameters(ar);

      int nAtomType = simulation().nAtomType();
      if (!prefactors_.isAllocated()) {
         prefactors_.allocate(nAtomType);
         cv_.allocate(nAtomType);
         cr_.allocate(nAtomType);
      }
      //  Note: Values of prefactors_, cv_, cr_ calculated in setup()
   }

   /*
   * Save internal state to an archive.
   */
   void NvtLangevinIntegrator::save(Serializable::OArchive &ar)
   {
      ar << dt_;
      ar << gamma_;
      Integrator::save(ar);
   }
 
   /*
   * Setup at beginning of run, before entering main loop.
   */ 
   void NvtLangevinIntegrator::setup()
   {
      // Precondition
      const EnergyEnsemble& energyEnsemble = simulation().energyEnsemble();
      if (!energyEnsemble.isIsothermal()) {
         UTIL_THROW("Energy ensemble is not isothermal");
      }

      // Initialize state and clear statistics on first usage.
      if (!isSetup()) {
         clear();
         setIsSetup();
      }

      // Exchange atoms, build pair list, compute forces.
      setupAtoms();

      // Set constants that are independent of atom type
      double cv = (exp(-dt_*gamma_) - 1.0)/dt_;
      double temp = energyEnsemble.temperature();
      double d = 2.0/(1.0 + exp(-dt_*gamma_));
      double cr = 12.0*temp*d*(1.0 - exp(-2.0*dt_*gamma_))/(dt_*dt_);

      // Loop over atom types
      double dtHalf = 0.5*dt_;
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = dtHalf/mass;
         cv_[i] = mass*cv;
         cr_[i] = sqrt(mass*cr);
      }

   }

   /*
   * First half of velocity-Verlet update.
   */
   void NvtLangevinIntegrator::integrateStep1()
   {
      Vector dv;
      Vector dr;
      double prefactor; // = 0.5*dt/mass
      AtomIterator atomIter;

      // 1st half of velocity Verlet.
      atomStorage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         prefactor = prefactors_[atomIter->typeId()];

         dv.multiply(atomIter->force(), prefactor);
         atomIter->velocity() += dv;

         dr.multiply(atomIter->velocity(), dt_);
         atomIter->position() += dr;
      }
   }

   /*
   * Second half of velocity-Verlet update.
   */
   void NvtLangevinIntegrator::integrateStep2()
   {
      Util::Random& random = simulation().random();
      Vector dv;
      Vector df;
      double cr;
      AtomIterator atomIter;
      int typeId, j;

      // 2nd half of velocity Verlet
      atomStorage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         typeId = atomIter->typeId();

         // Add Langevin drag and random force to atomic force
         df.multiply(atomIter->velocity(), cv_[typeId]);
         cr = cr_[typeId];
         for (j=0; j < Dimension; ++j) {
            df[j] += (random.uniform() - 0.5)*cr;
         }
         atomIter->force() += df;

         // Update velocity (half step)
         dv.multiply(atomIter->force(), prefactors_[typeId]);
         atomIter->velocity() += dv;
      }

      // Notify observers of change in velocity
      simulation().velocitySignal().notify();
   }

}
