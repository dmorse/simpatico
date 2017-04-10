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
#include <util/space/Vector.h>           // member template parameter
#include <util/space/Tensor.h>           // member template parameter
#include <util/containers/Pair.h>        // member template parameter
#include <util/containers/GArray.h>      // member template
#include <util/misc/Setable.h>           // member template
#include <mcMd/chemistry/AtomType.h>     // member template parameter
#include <util/containers/Array.h>       // member class template
#include <util/boundary/Boundary.h>      // typedef

#include <complex>

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
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      //@}
      /// \name System energy and stress.
      //@{
      
      /**
      * Generate wavevectors for the current boundary and kCutoff.
      */
      virtual void makeWaves();

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
   
      //@}
      /// \name Accessors 
      //@{

      /**
      * Unset KSpace energy.
      */
      void unsetEnergy();

      // Return K-Space contributions
      virtual double kSpaceEnergy(); 
      virtual double rSpaceEnergy();
      virtual double energy();

      /**
      * Unset KSpace stress.
      */
      void unsetStress();

      virtual Tensor kSpaceStress();
      virtual Tensor rSpaceStress();
      virtual Tensor stress() ;
      
      //@}
   

      /**
      * Current number of wavevectors with |k| < kCutoff.
      */
      int nWave() const;

      EwaldRSpaceAccumulator& rSpaceAccumulator()   
      {  return rSpaceAccumulator_; }

      EwaldInteraction& ewaldInteraction()
      { return ewaldInteraction_; }

   private:

      // Ewald Interaction - core Ewald computations
      EwaldInteraction ewaldInteraction_;

      // Pointer to parent Simulation
      Simulation* simulationPtr_;

      // Pointer to parent System
      System* systemPtr_;

      // Pointer to boundary of associated System.
      Boundary* boundaryPtr_;

      // Pointer to array of atom types
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

      /// Unit Matrix (constant).
      Tensor unitTensor_;

      /// Prefactor for self-interaction correction.
      double selfPrefactor_;

      // R-space energy and stress contributions
      EwaldRSpaceAccumulator rSpaceAccumulator_;

      // K-space part of Coulomb energy
      Setable<double> kSpaceEnergy_;

      // K-space part of Coulomb stress.
      Setable<Tensor> kSpaceStress_;

      /**
      * Calculate Fourier coefficients of charge density.
      */
      void computeKSpaceCharge();

   };

   inline void MdEwaldPotential::unsetEnergy() 
   {  kSpaceEnergy_.unset(); }

   inline double MdEwaldPotential::kSpaceEnergy()
   {  return kSpaceEnergy_.value(); }

   inline double MdEwaldPotential::rSpaceEnergy() 
   {  return rSpaceAccumulator_.rSpaceEnergy(); }

   inline double MdEwaldPotential::energy() 
   {  
      computeEnergy();
      return kSpaceEnergy_.value() + rSpaceAccumulator_.rSpaceEnergy();
   }

   inline void MdEwaldPotential::unsetStress()
   {  kSpaceStress_.unset(); }

   inline Tensor MdEwaldPotential::kSpaceStress()
   {  return kSpaceStress_.value(); }

   inline Tensor MdEwaldPotential::rSpaceStress()
   {  return rSpaceAccumulator_.rSpaceStress(); }

   inline Tensor MdEwaldPotential::stress() 
   { 
      Tensor temp = kSpaceStress_.value();
             temp += rSpaceAccumulator_.rSpaceStress();
      return temp;
   } 
      


} 
#endif

