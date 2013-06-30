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
#include <util/space/Vector.h>
#include <util/space/Tensor.h>


namespace McMd
{

   using namespace Util;


   /*
   * Constructor .
   */
   KSpaceCoulombPotential::KSpaceCoulombPotential(System& system)
    : systemPtr_(&system),
      boundaryPtr_(&system.boundary())
   {  setClassName("KSpaceCoulombPotential"); }

   /**
   * Destructor (does nothing)
   */
   KSpaceCoulombPotential::~KSpaceCoulombPotential()
   {}

   /// \name Initialization
   //@{

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
      gamma_.allocate(dimensions_);
   }

   /*
   * Compute values of k-squared and damping function
   */
   void KSpaceCoulombPotential::computeKSq()
   {}

   /*
   * Calculate the Fourier space contribution to the Coulomb energy.
   */
   double energy()
   {
      double total = 0.0;
      // .....
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
