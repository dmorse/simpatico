#ifndef MCMD_EWALD_COULOMB_POTENTIAL_CPP
#define MCMD_EWALD_COULOMB_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "EwaldCoulombPotential.h" 
//#include "EwaldCoulombPair.h" 
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
//#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/math/Constants.h>
#include <stdlib.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   EwaldCoulombPotential::EwaldCoulombPotential(System& system)
    : CoulombPotential(system),
      kCutoff_(1.0)
   {
      // Note: Don't setClassName - using "CoulombPotential" base class name
   }

   /*
   * Destructor (does nothing)
   */
   EwaldCoulombPotential::~EwaldCoulombPotential()
   {}

   /*
   * Read parameters and initialize.
   */
   void EwaldCoulombPotential::readParameters(std::istream& in)
   {
      CoulombPotential::readParameters(in);
      read<double>(in, "kCutoff", kCutoff_);
      isInitialized_ = true;
   }

   /*
   * Get cutfoff wavenumber for long range interaction.
   */
   inline int EwaldCoulombPotential::nWave() const
   {  return waves_.size(); }

   /*
   * Generate waves using kCutoff; allocate memories for associated variables.
   *
   * Comments:
   *   Only half of waves are stored.
   *   The first index is always non-negative.
   *
   */
   void EwaldCoulombPotential::makeWaves()
   {
      Vector    b0, b1, b2;    // Recprocal basis vectors.
      Vector    q0, q1, q;     // Partial and complete wavevectors.
      Vector    kv;            // Integer wavevector expressed as real vector
      double    prefactor(-0.25/alpha_/alpha_);
      double    kCutoffSq(kCutoff_*kCutoff_), ksq;
      double    pi2(2.0*Constants::Pi);
      IntVector maxK, k;       // Max and running wave indices.
      int       mink1, mink2;  // Minimum k-indices
      int       j;

      b0 = boundaryPtr_->reciprocalBasisVector(0);
      b1 = boundaryPtr_->reciprocalBasisVector(1);
      b2 = boundaryPtr_->reciprocalBasisVector(2);

      // Get max wave indices and reserve arrays.
      for (j=0; j < Dimension; ++j) {
         maxK[j] = ceil(kCutoff_*boundaryPtr_->bravaisBasisVector(j).abs()/pi2);
      }

      if (waves_.capacity() == 0) {
         int capacity; 
         capacity = ((2*maxK[0] + 1) * (2*maxK[1] + 1) * (2*maxK[2] + 1) - 1)/2;
         waves_.reserve(capacity);
         ksq_.reserve(capacity);
         g_.reserve(capacity);
         rho_.reserve(capacity);
      } else {
         waves_.clear();
         ksq_.clear();
         g_.clear();
         rho_.clear();
      }

      // Accumulate waves, and wave-related properties.
      q0.multiply(b0, -1);
      for (k[0] = 0; k[0] <= maxK[0]; ++k[0]) { // First index always non-negative.
         q0 += b0;

         mink1 = (k[0] == 0 ? 0 : -maxK[1]);
         q1.multiply(b1, mink1 - 1);
         q1 += q0;
         for (k[1] = mink1; k[1] <= maxK[1]; ++k[1]) {
            q1 += b1;

            mink2 = (k[0] == 0 && k[1] == 0 ? 1 : -maxK[2]);
            q.multiply(b2, mink2 - 1);
            q += q1;

            for (k[2] = mink2; k[2] <= maxK[2]; ++k[2]) {
               q += b2;

               ksq = q.square();
               if (ksq <= kCutoffSq) {
                  for (j = 0; j < Dimension; ++j) {
                    kv[j] = (double)k[j];
                  }
                  waves_.append(kv);
                  ksq_.append(ksq);
                  g_.append(exp(prefactor*ksq)/ksq);
               }

            } // for k[2]
         } // for k[1]
      } // for k[0]

      // Allocate fourier modes for charge density.
      rho_.resize(waves_.size());
   }

   /*
   * Calculate fourier modes of charge density.
   */
   void EwaldCoulombPotential::computeKSpaceCharge()
   {
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      Vector  rg;     // scaled atom position vector
      double  dotqr;  // dot product between q and r
      double  x, y;   // real and imaginary parts of phasor
      double  charge; // atom charge
      double  TwoPi;  // 2.0*pi
      int  i;         // array index

      TwoPi = 2.0*Constants::Pi;

      // Clear rho for all waves
      for (i = 0; i < rho_.size(); ++i) {
         rho_[i] = std::complex<double>(0.0, 0.0);
      }

      // Loop over species, molecules atoms
      int  nSpecies = simulationPtr_->nSpecies();
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               boundaryPtr_->transformCartToGen(atomIter->position(), rg);
               charge = (*atomTypesPtr_)[atomIter->typeId()].charge();

               // Loop over waves
               for (i = 0; i < waves_.size(); ++i) {
                  dotqr = rg.dot(waves_[i])*TwoPi;
                  x = charge*cos(dotqr);
                  y = charge*sin(dotqr);
                  rho_[i] += std::complex<double>(x, y);
               }

            } // For atoms.
         } // For molecules.
      } // For species.

   }

   /*
   * Calculate the k-space contribution to the Coulomb energy.
   */
   double EwaldCoulombPotential::kspaceEnergy()
   {
      double total = 0.0;
      double x, y;

      // Compute Fourier components of charge density.
      computeKSpaceCharge();

      for (int i = 0; i < waves_.size(); ++i) {
         x = rho_[i].real();
         y = rho_[i].imag();
         total += (x*x + y*y)*g_[i];
      }
      total *= 0.5 / (epsilon_*boundaryPtr_->volume());

      // Correct for conjugate wave contribution and return.
      return (2.0 * total);
   }

   /*
   * Add k-space Coulomb forces for all atoms.
   */
   void EwaldCoulombPotential::addKSpaceForces()
   {
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      Vector  b[Dimension];  // array of reciprocal basis vectors
      Vector  rg;            // scaled atom position vector
      Vector  fg;            // generalized force on atom
      Vector  df;            // force contribution
      double  dotqr;         // Dot product between q and r
      double  charge;        // atom charge
      double  x, y;          // Real and imaginary parts of phasor
      double  prefactor;     // -2/(epsilon*volume)
      double  TwoPi;         // 2.0*pi
      int  nSpecies(simulationPtr_->nSpecies());
      int  type;
      int  i, j;

      // Constants
      TwoPi   = 2.0*Constants::Pi;
      prefactor = -2.0/(epsilon_*boundaryPtr_->volume());

      // Compute Fourier components of charge density.
      computeKSpaceCharge();
    
      // Store reciprocal lattice vectors 
      for (j = 0; j < Dimension; ++j) {
         b[j] = boundaryPtr_->reciprocalBasisVector(j);
      }

      // Loop over species, molecules atoms
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               boundaryPtr_->transformCartToGen(atomIter->position(), rg);
               type = atomIter->typeId();
               charge = (*atomTypesPtr_)[type].charge();

               // Loop over waves
               fg.zero();
               for (i = 0; i < waves_.size(); ++i) {
                  df = waves_[i];
                  dotqr = rg.dot(df)*TwoPi;
                  x = cos(dotqr);
                  y = sin(dotqr);
                  df *= g_[i]*( x*rho_[i].imag() - y*rho_[i].real() );
                  fg += df;
               }
               fg *= charge*prefactor;

               // Transform to Cartesian coordinates
               for (j = 0; j < Dimension; ++j) {
                  df.multiply(b[j], fg[j]);
                  atomIter->force() += df;
               }

            } // for atoms
         } // for molecules
      } // for species

   }

   /*
   * Compute total nonbonded pressure
   */
   void EwaldCoulombPotential::computeKSpaceStress(double& stress)
   {}

   /*
   * Compute x, y, z nonbonded pressures.
   */
   void EwaldCoulombPotential::computeKSpaceStress(Util::Vector& stress)
   {}

   /*
   * Compute stress tensor.
   */
   void EwaldCoulombPotential::computeKSpaceStress(Util::Tensor& stress)
   {}
 
} 
#endif
