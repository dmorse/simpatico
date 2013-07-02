#ifndef MCMD_COULOMB_POTENTIAL_CPP
#define MCMD_COULOMB_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CoulombPotential.h" 
#include "EwaldCoulombPair.h" 
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
   CoulombPotential::CoulombPotential(System& system)
    : simulationPtr_(&system.simulation()),
      systemPtr_(&system),
      boundaryPtr_(&system.boundary()),
      pairInteractionPtr_(0),
      epsilon_(1.0),
      alpha_(1.0),
      rCutoff_(1.0),
      maxK_()
   {  setClassName("CoulombPotential"); }

   /*
   * Destructor (does nothing)
   */
   CoulombPotential::~CoulombPotential()
   {}

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
      read<double>(in, "alpha", alpha_);
      read<double>(in, "rCtuoff", rCutoff_);
      read<IntVector>(in, "maxK", maxK_);
      for (int i=0; i < Dimension; ++i) {
         dimensions_[i] = 2*maxK_[i] + 1;
      }
      rho_.allocate(dimensions_);
      ksq_.allocate(dimensions_);
      g_.allocate(dimensions_);
      pairInteractionPtr_->set(epsilon_, alpha_, rCutoff_);
   }

   /*
   * Compute values of k-squared and damping function
   */
   void CoulombPotential::computeKSq()
   {
      IntVector p;            // 0,...,dimensions_[i] - 1
      Vector    b0, b1, b2;   // recprocal basis vectors
      Vector    q0, q1, q;    // partial and complete wavevectors
      double    prefactor = -0.25/alpha_/alpha_;
      double    ksq;
      int       rank, zeroRank;

      b0 = boundaryPtr_->reciprocalBasisVector(0);
      b1 = boundaryPtr_->reciprocalBasisVector(1);
      b2 = boundaryPtr_->reciprocalBasisVector(2);

      // Calculate array rank for zero wavevector
      p[0] = maxK_[0];
      p[1] = maxK_[1];
      p[2] = maxK_[2];
      zeroRank = ksq_.rank(p);

      rank = -1;
      q0.multiply(b0, -maxK_[0]-1);
      for (p[0] = 0; p[0] < dimensions_[0]; ++p[0]) {
         q0 += b0;

         q1.multiply(b1, -maxK_[1]-1);
         q1 += q0;
         for (p[1] = 0; p[1] < dimensions_[1]; ++p[1]) {
            q1 += b1;

            q.multiply(b2, -maxK_[2]-1);
            q += q1;
            for (p[2] = 0; p[2] < dimensions_[2]; ++p[2]) {
               q += b2;
               ++rank;
               if (rank == zeroRank) {
                  ksq_(p) = 0.0;;
                  g_(p) = 0.0;
               } else {
                  ksq = q.square();
                  ksq_(p) = ksq;
                  g_(p) = exp(prefactor*ksq)/ksq;
               }
            } // for p[2]

         } // for p[1]

      } //  for p[0] 
   }

   /*
   * Calculate the k-space contribution to the Coulomb energy.
   */
   double CoulombPotential::kspaceEnergy()
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
