#ifndef MCMD_MD_EWALD_POTENTIAL_H
#define MCMD_MD_EWALD_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/coulomb/MdCoulombPotential.h>      // base class
#include <mcMd/potentials/coulomb/EwaldRSpaceAccumulator.h>  // member
#include <mcMd/potentials/coulomb/EwaldInteraction.h>        // member

#include <util/space/IntVector.h>        // member template parameter
#include <util/space/Tensor.h>           // member template parameter
#include <util/containers/Pair.h>        // member template parameter
#include <util/containers/GArray.h>      // member template
#include <util/misc/Setable.h>           // member template
#include <mcMd/chemistry/AtomType.h>     // member template parameter
#include <util/containers/Array.h>       // member class template
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

   typedef std::complex<double> DCMPLX;

   using namespace Util;

   /**
   * Ewald Coulomb potential class for MD simulations.
   *
   * This class implements the k-space sums for an Ewald  
   * summation of the Coulomb energy and forces.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class MdEwaldPotential : public MdCoulombPotential
   {

   public:

      /**
      * Constructor.
      */
      MdEwaldPotential(System& system);

      /**
      * Destructor (does nothing).
      */
      virtual ~MdEwaldPotential();

      /// \name Initialization
      //@{

      /**
      * Read parameters and initialize.
      *
      * \param in input stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Generate wavevectors for the current boundary and kCutoff.
      */
      virtual void makeWaves();

      //@}
      /// \name System energy and stress.
      //@{
      
      /**
      * Add k-space Coulomb forces for all atoms.
      */
      virtual void addForces();

      /**
      * Calculate the long range kspace part of Coulomb energy.
      */
      virtual void computeEnergy();

      /**
      * Compute kspace part of Coulomb pressure.
      *
      * \param stress (output) pressure
      */
      virtual void computeStress();
   
      // Unset KSpace energy and stress components
      virtual void unsetEnergy();
      virtual void unsetStress();

      //@}
      /// \name Accessors (const)
      //@{

      // Return K-Space contributions
      virtual double kSpaceEnergy() const;
      //virtual double kSpacePressure() const;
      virtual Tensor kSpaceStress() const;

      // Return R-Space contributions
      virtual double rSpaceEnergy() const;
      //virtual double rSPacePressure() const;
      virtual Tensor rSpaceStress() const;

      // Return total Coulomb energy and stress (kspace + rspace)
      virtual double energy() ;
      //virtual double pressure() ;
      virtual Tensor stress() ;
      

      /**
      * Current number of wavevectors with |k| < kCutoff.
      */
      int nWave() const;

      //@}
      
      // R-space energy and stress contributions
      EwaldRSpaceAccumulator rSpaceAccumulator_;

      // Ewald Interaction - core Ewald computations
      EwaldInteraction ewaldInteraction_;

   protected:

      Simulation* simulationPtr_;

      System* systemPtr_;

      Boundary* boundaryPtr_;

      const Array<AtomType>* atomTypesPtr_;

      /// Exponential factor accessors.
      double  base0_, base1_, base2_;
      double  upper0_, upper1_, upper2_;
      GArray<DCMPLX> fexp0_;
      GArray<DCMPLX> fexp1_;
      GArray<DCMPLX> fexp2_;

      /// Wave vector indices.
      GArray<IntVector> waves_;

      //real space vector indices.
      GArray<Vector> reals_;

      /// Values of square of Fourier wavevector.
      GArray<double> ksq_;

      /// Regularized Green's function (Gaussian/ksq)
      GArray<double> g_;

      /// Fourier modes of charge density.
      GArray<DCMPLX> rho_;

      /**
      * Calculate fourier modes of charge density.
      */
      void computeKSpaceCharge();

   private:

      // K-space part of Coulomb energy
      Setable<double> kSpaceEnergy_;

      // K-space part of Coulomb stress.
      Setable<Tensor> kSpaceStress_;

      /// Prefactor for self-interaction correction.
      double selfPrefactor_;

      /// Unit Matrix (constant).
      Tensor unitTensor_;

   };

   inline void MdEwaldPotential::unsetEnergy() 
   {  kSpaceEnergy_.unset(); }

   inline void MdEwaldPotential::unsetStress()
   {  kSpaceStress_.unset(); }

   inline double MdEwaldPotential::kSpaceEnergy() const
   {  return kSpaceEnergy_.value(); }

   inline Tensor MdEwaldPotential::kSpaceStress() const
   {  return kSpaceStress_.value(); }

   inline double MdEwaldPotential::rSpaceEnergy() const
   {  return rSpaceAccumulator_.rSpaceEnergy(); }

   inline Tensor MdEwaldPotential::rSpaceStress() const
   {  return rSpaceAccumulator_.rSpaceStress(); }

   inline double MdEwaldPotential::energy() 
   {  computeEnergy();
      return kSpaceEnergy_.value() + rSpaceAccumulator_.rSpaceEnergy();
   }

   //inline double pressure()
   //{  return kSpace}

   inline Tensor MdEwaldPotential::stress() 
   { 
      Tensor temp = kSpaceStress_.value();
             temp += rSpaceAccumulator_.rSpaceStress();
      return temp;
   } 
      


} 
#endif


