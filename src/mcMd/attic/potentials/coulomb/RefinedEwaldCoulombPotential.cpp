#ifndef MCMD_REFINED_EWALD_COULOMB_POTENTIAL_CPP
#define MCMD_REFINED_EWALD_COULOMB_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "RefinedEwaldCoulombPotential.h" 
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/math/Constants.h>
#include <stdlib.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   RefinedEwaldCoulombPotential::RefinedEwaldCoulombPotential(System& system)
    : CoulombPotential(),
      CoulombSystemMixIn(system),
      kCutoff_(1.0)
   {
      // Note: Don't setClassName - using "CoulombPotential" base class name
   }

   /*
   * Destructor (does nothing)
   */
   RefinedEwaldCoulombPotential::~RefinedEwaldCoulombPotential()
   {}

   /*
   * Read parameters and initialize.
   */
   void RefinedEwaldCoulombPotential::readParameters(std::istream& in)
   {
      CoulombPotential::readParameters(in);
      read<double>(in, "kCutoff", kCutoff_);
      isInitialized_ = true;
   }

   /*
   * Get cutfoff wavenumber for long range interaction.
   */
   inline int RefinedEwaldCoulombPotential::nWave() const
   {  return waves_.size(); }

   /*
   * Generate waves using kCutoff; allocate memories for associated variables.
   *
   * Comments:
   *   Only half of waves are stored.
   *   The first index is always non-negative.
   *
   */
   void RefinedEwaldCoulombPotential::makeWaves()
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
      int       upper0, upper1, upper2; // Upper wave numbers.

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

         fexp0_.reserve(maxK[0] + 1);
         fexp1_.reserve(2*maxK[1] + 1);
         fexp2_.reserve(2*maxK[2] + 1);
      } else {
         waves_.clear();
         ksq_.clear();
         g_.clear();
         rho_.clear();

         fexp0_.clear();
         fexp1_.clear();
         fexp2_.clear();
      }

      // Accumulate waves, and wave-related properties.
      base0_ = 0;
      upper0 = -maxK[0];

      base1_ = maxK[1];
      upper1 = -base1_;

      base2_ = maxK[2];
      upper2 = -base2_;

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

               ksq = double(q.square());
               if (ksq <= kCutoffSq) {

                  if (k[0] > upper0) upper0 = k[0];

                  if (k[1] < base1_) base1_ = k[1];
                  if (k[1] > upper1) upper1 = k[1];

                  if (k[2] < base2_) base2_ = k[2];
                  if (k[2] > upper2) upper2 = k[2];

                  waves_.append(k);
                  ksq_.append(ksq);
                  g_.append(exp(prefactor*ksq)/ksq);
               }

            } // for k[2]
         } // for k[1]
      } // for k[0]

      fexp0_.resize(upper0 - base0_ + 1);
      fexp1_.resize(upper1 - base1_ + 1);
      fexp2_.resize(upper2 - base2_ + 1);

      // Allocate fourier modes for charge density.
      rho_.resize(waves_.size());
   }

   /*
   * Calculate fourier modes of charge density.
   */
   void RefinedEwaldCoulombPotential::computeKSpaceCharge()
   {
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      Vector  rg;     // scaled atom position vector
      IntVector q;    // wave indices
      Vector  qv;     // wave vector
      double  dotqr;  // dot product between q and r
      double  x, y;   // real and imaginary parts of phasor
      double  charge; // atom charge
      double  TwoPi;  // 2.0*pi
      int  i, j;      // array index

      TwoPi = 2.0*Constants::Pi;

      // Clear rho for all waves
      for (i = 0; i < rho_.size(); ++i) {
         rho_[i] = DCMPLX(0.0, 0.0);
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
                  q = waves_[i];
                  for (j = 0; j < Dimension; ++j) {
                     qv[j] = double(q[j]);
                  }
                  dotqr = rg.dot(qv)*TwoPi;
                  x = charge*cos(dotqr);
                  y = charge*sin(dotqr);
                  rho_[i] += DCMPLX(x, y);
               }

            } // For atoms.
         } // For molecules.
      } // For species.

   }

   /*
   * Calculate the k-space contribution to the Coulomb energy.
   */
   double RefinedEwaldCoulombPotential::kspaceEnergy()
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
   void RefinedEwaldCoulombPotential::addKSpaceForces()
   {
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      Vector  b[Dimension];  // array of reciprocal basis vectors
      Vector  rg;            // scaled atom position vector
      Vector  fg;            // generalized force on atom
      Vector  df;            // force contribution
      IntVector q;           // wave index
      double  charge;        // atom charge
      double  x, y;          // Real and imaginary parts of phasor
      double  prefactor;     // -2/(epsilon*volume)
      double  TwoPi;         // 2.0*pi
      DCMPLX  TwoPiIm;       // 2.0*pi*I
      double  EPS(1.0E-10);  // Tiny number to check if is charge
      int  nSpecies(simulationPtr_->nSpecies());
      int  type;
      int  i, j;
      DCMPLX  de, expfactor;

      // Constants
      TwoPi   = 2.0*Constants::Pi;
      TwoPiIm = TwoPi * Constants::Im;
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
               type = atomIter->typeId();
               charge = (*atomTypesPtr_)[type].charge();

               if (fabs(charge) > EPS) {
                  boundaryPtr_->transformCartToGen(atomIter->position(), rg);

                  // Tabulate the exponential factors.
                  fexp0_[0] = exp(TwoPiIm * rg[0] * double(base0_));
                  de = exp(TwoPiIm * rg[0]);
                  for (i = 1; i < fexp0_.size(); ++i) {
                     fexp0_[i] = fexp0_[i-1] * de;
                  }

                  fexp1_[0] = exp(TwoPiIm * rg[1] * double(base1_));
                  de = exp(TwoPiIm * rg[1]);
                  for (i = 1; i < fexp1_.size(); ++i) {
                     fexp1_[i] = fexp1_[i-1] * de;
                  }

                  fexp2_[0] = exp(TwoPiIm * rg[2] * double(base2_));
                  de = exp(TwoPiIm * rg[2]);
                  for (i = 1; i < fexp2_.size(); ++i) {
                     fexp2_[i] = fexp2_[i-1] * de;
                  }

                  // Accumulating forces.
                  fg.zero();
                  for (i = 0; i < waves_.size(); ++i) {
                     q = waves_[i];
                     for (j = 0; j < Dimension; ++j) {
                        df[j] = double(q[j]);
                     }
                     expfactor = fexp0_[q[0]-base0_] * fexp1_[q[1]-base1_] * fexp2_[q[2]-base2_];
                     x = expfactor.real();
                     y = expfactor.imag();
                     df *= g_[i]*( x*rho_[i].imag() - y*rho_[i].real() );
                     fg += df;
                  }
                  fg *= charge*prefactor;

                  // Transform to Cartesian coordinates
                  for (j = 0; j < Dimension; ++j) {
                     df.multiply(b[j], fg[j]);
                     atomIter->force() += df;
                  }
               }

            } // for atoms
         } // for molecules
      } // for species

   }

   /*
   * Compute total nonbonded pressure
   */
   void RefinedEwaldCoulombPotential::computeKSpaceStress(double& stress)
   {}

   /*
   * Compute x, y, z nonbonded pressures.
   */
   void RefinedEwaldCoulombPotential::computeKSpaceStress(Util::Vector& stress)
   {}

   /*
   * Compute stress tensor.
   */
   void RefinedEwaldCoulombPotential::computeKSpaceStress(Util::Tensor& stress)
   {}
 
} 
#endif
