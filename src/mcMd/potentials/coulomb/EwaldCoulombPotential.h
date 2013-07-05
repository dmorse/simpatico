#ifndef MCMD_EWALD_COULOMB_POTENTIAL_H
#define MCMD_EWALD_COULOMB_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CoulombPotential.h"            // base class
#include <mcMd/chemistry/AtomType.h>     // Member template parameter
#include <util/space/Vector.h>           // member template parameter
#include <util/containers/GArray.h>      // member template
#include <util/boundary/Boundary.h>      // typedef

#include <complex>

namespace Util
{
   class Vector;
   class Tensor;
}

namespace McMd
{

   class Simulation;
   class System;
   class EwaldCoulombPair;

   using namespace Util;

   /**
   * Ewald Coulomb potential class.
   *
   * This class implements the k-space sums for an Ewald  
   * summation of the Coulomb forces.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class EwaldCoulombPotential : public CoulombPotential
   {

   public:

      /**
      * Constructor.
      */
      EwaldCoulombPotential(System& system);

      /**
      * Destructor (does nothing).
      */
      virtual ~EwaldCoulombPotential();

      /// \name Initialization
      //@{

      /**
      * Read parameters and initialize.
      *
      * \param in input stream
      */
      void readParameters(std::istream& in);

      /**
      * Generate wavevectors for the current boundary and kCutoff.
      */
      virtual void makeWaves();

      //@}
      /// \name System energy and stress.
      //@{

      /**
      * Calculate the long range kspace part of Coulomb energy.
      */
      virtual double kspaceEnergy();

      /**
      * Add k-space Coulomb forces for all atoms.
      */
      virtual void addKSpaceForces();

      /**
      * Compute kspace part of Coulomb pressure.
      *
      * \param stress (output) pressure
      */
      virtual void computeKSpaceStress(double& stress);

      /**
      * Compute kspace part of Coulomb x, y, z pressures.
      *
      * \param stress (output) pressures
      */
      virtual void computeKSpaceStress(Util::Vector& stress);

      /**
      * Compute kspace part of Coulomb stress tensor.
      *
      * \param stress (output) stress tensor
      */
      virtual void computeKSpaceStress(Util::Tensor& stress);
    
      //@}
      /// \name Accessors (const)
      //@{

      /**
      * Current number of wavevectors with |k| < kCutoff.
      */
      int nWave() const;

      //@}

   private:

      /// Cutoff wavenumber for long range interaction.
      double kCutoff_;

      /// Wave vector indices.
      GArray<Vector> waves_;

      /// Values of square of Fourier wavevector.
      GArray<double> ksq_;

      /// Regularized Green's function (Gaussian/ksq)
      GArray<double> g_;

      /// Fourier modes of charge density.
      GArray< std::complex<double> > rho_;

      /**
      * Calculate fourier modes of charge density.
      */
      void computeChargeKMode();

   };

} 
#endif
