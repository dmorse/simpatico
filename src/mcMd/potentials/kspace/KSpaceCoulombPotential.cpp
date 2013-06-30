#ifndef MCMD_KSPACE_COULOMB_POTENTIAL_CPP
#define MCMD_KSPACE_COULOMB_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "KSpaceCoulombPotential.h" 
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor .
   */
   KSpaceCoulombPotential::KSpaceCoulombPotential(System& system)
    : 
      simulationPtr_(&system.simulation()),
      systemPtr_(&system),
      boundaryPtr_(&system.boundary())
   {  setClassName("KSpaceCoulombPotential"); }

   /*
   * Destructor (does nothing)
   */
   KSpaceCoulombPotential::~KSpaceCoulombPotential()
   {}

   /*
   * Read parameters and initialize.
   */
   void KSpaceCoulombPotential::readParam(std::istream& in)
   {
      read<double>(in, "epsilon", epsilon_);
      read<double>(in, "alpha", alpha_);
      read<IntVector>(in, "maxL", maxL_);
      for (int i=0; i < Dimension; ++i) {
         dimensions_[i] = 2*maxL_[i] + 1;
      }
      rho_.allocate(dimensions_);
      ksq_.allocate(dimensions_);
      g_.allocate(dimensions_);
   }

   /*
   * Compute values of k-squared and damping function
   */
   void KSpaceCoulombPotential::computeKSq()
   {
      IntVector p;            // 0,...,dimensions_[i] - 1
      IntVector k;            // -maxL_[i], ...., maxL_[i]
      Vector    b0, b1, b2;   // recprocal basis vectors
      Vector    q0, q1, q;   // partial and complete wavevectors
      double    prefactor = -0.25/alpha_;
      double    ksq;
      int       rank, zeroRank;

      b0 = boundaryPtr_->reciprocalBasisVector(0);
      b1 = boundaryPtr_->reciprocalBasisVector(1);
      b2 = boundaryPtr_->reciprocalBasisVector(2);

      // Calculate array rank for zero wavevector
      p[0] = maxL_[0];
      p[1] = maxL_[1];
      p[2] = maxL_[2];
      zeroRank = ksq_.rank(p);

      rank = 0;
      k[0] = -maxL_[0];
      for (p[0] = 0; p[0] < dimensions_[0]; ++p[0], k[0]) {
         q0.multiply(b0, k[0]);

         k[1] = -maxL_[1];
         for (p[1] = 0; p[1] < dimensions_[1]; ++p[1], ++k[1]) {
            q1.multiply(b1, k[1]);
            q1 += q0;

            q.multiply(b2, -maxL_[2]);
            q += q1;
            k[2] = -maxL_[2];
            for (p[2] = 0; p[2] < dimensions_[2]; ++p[2], k[2]) {
               if (rank == zeroRank) {
                  ksq_(p) = 0.0;;
                  g_(p) = 0.0;
               } else {
                  ksq = q.square();
                  ksq_(p) = ksq;
                  g_(p) = exp(prefactor*ksq)/ksq;
               }
               q += b2;
               ++rank;
            } // for p[2]

         } // for p[1]

      } //  for p[0] 
   }

   /*
   * Calculate the Fourier space contribution to the Coulomb energy.
   */
   double KSpaceCoulombPotential::energy()
   {
      Vector r;
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      int nSpecies = simulationPtr_->nSpecies();

      double total = 0.0;
      for (int iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               // Do something
            }
         }
      }
      return total;
   }

   #if 0
   /*
   * Compute total nonbonded pressure
   */
   void KSpaceCoulombPotential::computeStress(double& stress) const
   {}

   /*
   * Compute x, y, z nonbonded pressures.
   */
   void KSpaceCoulombPotential::computeStress(Util::Vector& stress) const 
   {}

   /*
   * Compute stress tensor.
   */
   void KSpaceCoulombPotential::computeStress(Util::Tensor& stress) const
   {}
   #endif
 
} 
#endif
