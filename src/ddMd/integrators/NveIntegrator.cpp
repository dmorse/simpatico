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
    : TwoStepIntegrator(simulation)
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
      #if 0
      atomStorage().clearSnapshot();
      exchanger().exchange();
      pairPotential().buildCellList();
      if (!UTIL_ORTHOGONAL) {
         atomStorage().transformGenToCart(boundary());
      }
      atomStorage().makeSnapshot();
      pairPotential().buildPairList();
      simulation().computeForces();
      #endif
      // Exchange atoms, build pair list, compute forces.
      setupAtoms();

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

   void NveIntegrator::integrateStep1()
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

   void NveIntegrator::integrateStep2()
   {
      Vector dv;
      double prefactor; // = 0.5*dt/mass
      AtomIterator atomIter;

      // 2nd half of velocity Verlet
      atomStorage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         prefactor = prefactors_[atomIter->typeId()];
         dv.multiply(atomIter->force(), prefactor);
         atomIter->velocity() += dv;
      }

   }

}
#endif
