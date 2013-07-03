#ifndef MCMD_MD_COULOMB_POTENTIAL_CPP
#define MCMD_MD_COULOMB_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include "MdCoulombPotential.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MdCoulombPotential::MdCoulombPotential(System& system)
      : CoulombPotential(system)
   { setClassName("MdCoulombPotential"); }

   /**
   * Destructor (does nothing)
   */
   MdCoulombPotential::~MdCoulombPotential()
   {}

   /*
   * Add k-space Coulomb forces for all atoms.
   */
   void MdCoulombPotential::addKSpaceForces()
   {
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      int       nSpecies(simulationPtr_->nSpecies());
      Vector    rg;            // General atom position vector.
      Vector    b0, b1, b2;    // Recprocal basis vectors.
      Vector    vq, vqtmp;     // Wavevector.
      Vector    force;         // Force.
      IntVector q;             // Wavenumber.
      double    dotqr;         // Dot product between q and r.
      double    x, y, fcoeff;

      b0 = boundaryPtr_->reciprocalBasisVector(0);
      b1 = boundaryPtr_->reciprocalBasisVector(1);
      b2 = boundaryPtr_->reciprocalBasisVector(2);

      // Comput Fourier components of charge density.
      computeChargeKMode();

      for (int i = 0; i < waves_.size(); ++i) {

         q = waves_[i];

         vq.multiply(b0, q[0]);
         vqtmp.multiply(b1, q[1]);
         vq += vqtmp;
         vqtmp.multiply(b2, q[2]);
         vq += vqtmp;
         vq *= -2.0;

         x = rho_[i].real();
         y = rho_[i].imag();

         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            systemPtr_->begin(iSpecies, molIter); 
            for ( ; molIter.notEnd(); ++molIter) {
               for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {

                  boundaryPtr_->transformCartToGen(atomIter->position(), rg);
                  dotqr = rg[0]*q[0] + rg[1]*q[1] + rg[2]*q[2];
                  fcoeff = y * cos(dotqr) - x * sin(dotqr);

                  force.multiply(vq, fcoeff);
                  atomIter->force() += force;

               } // For atoms.
            } // For molecules.
         } // For species.

      }
   }

} 
#endif
