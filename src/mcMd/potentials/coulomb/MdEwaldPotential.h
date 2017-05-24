#ifndef MCMD_MD_EWALD_POTENTIAL_H
#define MCMD_MD_EWALD_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/coulomb/MdCoulombPotential.h>       // base class
#include <mcMd/potentials/coulomb/EwaldRSpaceAccumulator.h>   // member
#include <simp/interaction/coulomb/EwaldInteraction.h>        // member

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
   using namespace Simp;

   /**
   * Ewald Coulomb potential class for MD simulations.
   *
   * This class implements the k-space sums in the Ewald
   * method for computing the Coulomb energy and forces.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class MdEwaldPotential : public MdCoulombPotential
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system.
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
      /// \name Parameters (get/set)
      //@{

      /**
      * Set a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param value  new value of parameter
      */
      void set(std::string name, double value);

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      */
      double get(std::string name) const;

      //@}
      /// \name System energy and stress.
      //@{

      /**
      * Generate wavevectors for the current boundary and kCutoff.
      */
      virtual void makeWaves();

      /**
      * Current number of wavevectors with |k| < kCutoff.
      */
      int nWave() const;

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
      */
      virtual void computeStress();

      //@}
      /// \name Miscellaneous Accessors
      //@{

      EwaldRSpaceAccumulator& rSpaceAccumulator()
      {  return rSpaceAccumulator_; }

      EwaldInteraction& ewaldInteraction()
      { return ewaldInteraction_; }

      //@}

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
      GArray<IntVector> intWaves_;

      //real space vector indices.
      GArray<Vector> reals_;

      /// Values of square of Fourier wavevector.
      GArray<double> ksq_;

      /// Regularized Green's function (Gaussian/ksq)
      GArray<double> g_;

      /// Fourier modes of charge density.
      GArray<DCMPLX> rho_;

      /// cutoff distance in k space
      double kSpaceCutoff_;

      /*
      * Calculate Fourier coefficients of charge density.
      */
      void computeKSpaceCharge();

   };

}
#endif

