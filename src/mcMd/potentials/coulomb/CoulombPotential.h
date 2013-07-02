#ifndef MCMD_COULOMB_POTENTIAL_H
#define MCMD_COULOMB_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <util/space/IntVector.h>        // members
#include <util/containers/GridArray.h>   // members
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

   using namespace Util;

   /**
   * Coulomb potential base class.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class CoulombPotential : public ParamComposite
   {

   public:

      /**
      * Constructor .
      */
      CoulombPotential(System& system);

      /**
      * Destructor (does nothing)
      */
      virtual ~CoulombPotential();

      /// \name Initialization
      //@{

      /**
      * Read parameters and initialize.
      *
      * \param in input stream
      */
      void readParam(std::istream& in);

      /**
      * Compute values of k-squared and influence function
      */
      void computeKSq();

      //@}

      /// \name System energy and stress.
      //@{

      /**
      * Calculate the long range kspace part of Coulomb energy.
      */
      double kspaceEnergy();

      #if 0
      /**
      * Compute kspace part of Coulomb pressure
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
      #endif
    
      //@}

      /**
      * Get value of inverse range parameter.
      */
      double alpha() const;

      /**
      * Get value of permittivity parameter.
      */
      double epsilon() const;

   protected:

      /// Fourier modes of charge density
      GridArray< std::complex<double> > rho_;

      /// Values of square of Fourier wavevector
      GridArray<double> ksq_;

      /// Regularized Green's function (Gaussian/ksq)
      GridArray<double> g_;

      /// Pointer to associated Simulation
      Simulation* simulationPtr_;

      /// Pointer to associated System
      System* systemPtr_;

      /// Pointer to associated Boundary
      Boundary* boundaryPtr_;

      /// Inverse range parameter
      double alpha_;
      
      /// Dielectric constant / permittivity
      double epsilon_;

      /// Maximum Miller indices for Fourier modes
      IntVector maxK_;

      /// Dimensions of Fourier mode grid = 2*maxK_ + 1
      IntVector dimensions_;

      /// Has readParam been called?
      bool isInitialized_;

   };

   // Inline functions

   /*
   * Get value of inverse range parameter.
   */
   inline double CoulombPotential::alpha() const
   {  return alpha_; }

   /*
   * Get value of epsilon (dielectric) parameter.
   */
   inline double CoulombPotential::epsilon() const
   {  return epsilon_; }

} 
#endif
