#ifndef MCMD_COULOMB_POTENTIAL_CPP
#define MCMD_COULOMB_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <stdlib.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/math/Constants.h>
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include "CoulombPotential.h" 
#include "EwaldCoulombPair.h" 

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   CoulombPotential::CoulombPotential(System& system)
    : simulationPtr_(&system.simulation()),
      systemPtr_(&system),
      boundaryPtr_(&system.boundary()),
      pairInteractionPtr_(0),
      epsilon_(1.0),
      alpha_(1.0),
      rCutoff_(1.0),
      kCutoff_(1.0),
      atomTypesPtr_(NULL),
      isInitialized_(false)
   { setClassName("CoulombPotential"); }

   /*
   * Destructor (does nothing)
   */
   CoulombPotential::~CoulombPotential()
   {}

   /*
   * Set pointer to array of AtomTypes.
   */
   void CoulombPotential::setAtomTypes(const Array<AtomType>& atomTypes)
   {  atomTypesPtr_ = &atomTypes; }


   /*
   * Destructor (does nothing)
   */
   void CoulombPotential::setPairInteraction(EwaldCoulombPair& pairInteraction)
   {
      pairInteractionPtr_ = &pairInteraction;
      pairInteractionPtr_->set(epsilon_, alpha_, rCutoff_);
   }

   /*
   * Read parameters and initialize.
   */
   void CoulombPotential::readParam(std::istream& in)
   {
      read<double>(in, "epsilon", epsilon_);
      read<double>(in, "alpha",   alpha_);
      read<double>(in, "rCutoff", rCutoff_);
      read<double>(in, "kCutoff", kCutoff_);
      pairInteractionPtr_->set(epsilon_, alpha_, rCutoff_);

      makeWaves();
      isInitialized_ = true;
   }

   /*
   * Generate waves using kCutoff; allocate memories for associated variables.
   *
   * Comments:
   *   Only half of waves are stored.
   *   The first index is always non-negative.
   *
   */
   void CoulombPotential::makeWaves()
   {
      Vector    b0, b1, b2;    // Recprocal basis vectors.
      IntVector maxK, k;       // Wave indices.
      int       mink1, mink2;  // Minimum k-indices.
      Vector    q0, q1, q;     // Partial and complete wavevectors.
      double    prefactor(-0.25/alpha_/alpha_);
      double    kCutoffSq(kCutoff_*kCutoff_), ksq;

      b0 = boundaryPtr_->reciprocalBasisVector(0);
      b1 = boundaryPtr_->reciprocalBasisVector(1);
      b2 = boundaryPtr_->reciprocalBasisVector(2);

      // Max wave indices.
      maxK[0] = ceil(kCutoff_ * boundaryPtr_->length(0) / 2.0 / Constants::Pi);
      maxK[1] = ceil(kCutoff_ * boundaryPtr_->length(1) / 2.0 / Constants::Pi);
      maxK[2] = ceil(kCutoff_ * boundaryPtr_->length(2) / 2.0 / Constants::Pi);

      // Accumulate waves, and wave-related properties.
      q0.multiply(b0, -1);
      for (k[0] = 0; k[0] <= maxK[0]; ++k[0]) { // First index always non-negative.
         q0 += b0;

         q1.multiply(b1, -maxK[1]-1);
         q1 += q0;

         mink1 = (k[0] == 0 ? 0 : -maxK[1]);
         for (k[1] = mink1; k[1] <= maxK[1]; ++k[1]) {
            q1 += b1;

            q.multiply(b2, -maxK[2]-1);
            q += q1;

            mink2 = (k[0] == 0 && k[1] == 0 ? 0 : -maxK[2]);
            for (k[2] = mink2; k[2] <= maxK[2]; ++k[2]) {

               if (k[0] + abs(k[1]) + abs(k[2]) > 0) {
                  ksq = q.square();
                  if (ksq <= kCutoffSq) {
                     waves_.append(k);
                     ksq_.append(ksq);
                     g_.append(exp(prefactor*ksq)/ksq);
                  }
               } else {
                  waves_.append(k);
                  ksq_.append(0.0);
                  g_.append(0.0);
               }

            } // for k[2]
         } // for k[1]
      } // for k[0]

      // Allocate fourier modes for charge density.
      rho_.resize(waves_.size());
   }

   /*
   * Compute values of k-squared and damping function.
   * Needed when box is reshaped.
   */
   void CoulombPotential::computeKSq()
   {
      Vector    b0, b1, b2;   // recprocal basis vectors
      Vector    q, qtmp;      // partial and complete wavevectors
      double    prefactor = -0.25/alpha_/alpha_;
      double    ksq;

      b0 = boundaryPtr_->reciprocalBasisVector(0);
      b1 = boundaryPtr_->reciprocalBasisVector(1);
      b2 = boundaryPtr_->reciprocalBasisVector(2);

      for (int i = 0; i < waves_.size(); ++i) {
         q.multiply(b0, waves_[i][0]);
         qtmp.multiply(b1, waves_[i][1]);
         q += qtmp;
         qtmp.multiply(b2, waves_[i][2]);
         q += qtmp;

         ksq = q.square();
         ksq_[i] = ksq;
         g_[i] = exp(prefactor*ksq) / ksq;
      }
   }

   /*
   * Calculate fourier modes of charge density.
   */
   void CoulombPotential::computeChargeKMode()
   {
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      int     nSpecies(simulationPtr_->nSpecies());
      std::complex<double> img(0.0, 1.0); // Imaginary number unit.
      Vector     rg;         // General atom position vector.
      IntVector  q;          // Wavenumber.
      int        i;          // Index for waves.
      int        type;       // Atom type id.
      double     dotqr;      // Dot product between q and r.

      for (i = 0; i < rho_.size(); ++i)
         rho_[i] = std::complex<double>(0.0, 0.0);

      for (int i = 0; i < waves_.size(); ++i) {
         q = waves_[i];

         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            systemPtr_->begin(iSpecies, molIter); 
            for ( ; molIter.notEnd(); ++molIter) {
               for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {

                  boundaryPtr_->transformCartToGen(atomIter->position(), rg);
                  dotqr = rg[0]*q[0] + rg[1]*q[1] + rg[2]*q[2];
                  type = atomIter->typeId();
                  rho_[i] += (*atomTypesPtr_)[type].charge() * exp(img*dotqr);

               } // For atoms.
            } // For molecules.
         } // For species.

      }

   }


   /*
   * Calculate the k-space contribution to the Coulomb energy.
   */
   double CoulombPotential::kspaceEnergy()
   {
      double total(0.0);
      double x, y;

      // Comput Fourier components of charge density.
      computeChargeKMode();

      for (int i = 0; i < waves_.size(); ++i) {
         x = rho_[i].real();
         y = rho_[i].imag();
         total += (x*x + y*y)*g_[i];
      }
      total *= 0.5 / epsilon_ / boundaryPtr_->volume();

      // Correct for conjugate wave contribution and return.
      return (2.0 * total);
   }

   /*
   * Add k-space Coulomb forces for all atoms.
   */
   void CoulombPotential::addKSpaceForces()
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

   #if 0
   /*
   * Compute total nonbonded pressure
   */
   void CoulombPotential::computeKSpaceStress(double& stress) const
   {}

   /*
   * Compute x, y, z nonbonded pressures.
   */
   void CoulombPotential::computeKSpaceStress(Util::Vector& stress) const 
   {}

   /*
   * Compute stress tensor.
   */
   void CoulombPotential::computeKSpaceStress(Util::Tensor& stress) const
   {}
   #endif
 
} 
#endif
