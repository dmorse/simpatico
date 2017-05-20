#ifndef MD_COULOMB_POTENTIAL_H
#define MD_COULOMB_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>                   // base class
#include <mcMd/potentials/coulomb/EwaldRSpaceAccumulator.h>  // member
#include <util/misc/Setable.h>                           // member template
#include <util/space/Tensor.h>                           // template parameter

#include <complex>

namespace McMd
{

   using namespace Util;

   /**
   * Coulomb potential for an Md simulation.
   *
   * This class computes the long-range k-space part of the
   * Coulomb forces and energies in a molecular dynamics (MD)
   * simulation, and provides accessors for both r-space and 
   * k-space contributions to the Coulomb energy and stress.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class MdCoulombPotential : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      MdCoulombPotential();

      /**
      * Destructor (does nothing).
      */
      virtual ~MdCoulombPotential();

      /**
      * Modify an interaction parameter, identified by a string.
      *
      * \param name  parameter name
      * \param value new value of parameter
      */
      virtual void set(std::string name, double value);
   

      /**
      * Get an interaction parameter value, identified by a string.
      *
      * \param name parameter name
      */
      virtual double get(std::string name) const;

      /// \name Waves (data that depends on Boundary).
      //@{

      /**
      * Are wavevectors and k-space influence function up to date?
      */
      bool hasWaves();

      /**
      * Generate wavevectors and influence function for this boundary.
      */
      virtual void makeWaves() = 0;

      /**
      * Unset all data that depends on the Boundary.
      *
      * Unsets stored values of waves, influence function, energy and stress.
      */
      void unsetWaves();

      /**
      * Current number of wavevectors.
      */
      virtual int nWave() const = 0;

      //@}
      /// \name Forces and Energy
      //@{

      /**
      * Add k-space Coulomb forces to forces on all atoms.
      */
      virtual void addForces() = 0;

      /**
      * Calculate the long range kspace part of Coulomb energy.
      */
      virtual void computeEnergy() = 0;

      /**
      * Unset k-space energy.
      */
      virtual void unsetEnergy();

      /**
      * Get long-range k-space part of Coulomb energy.
      * 
      * Recomputes iff necessary (i.e., if not set).
      */
      double kSpaceEnergy();

      /**
      * Return short-range r-space part of Coulomb energy.
      *
      * Uses associated pair potential to recompute iff necessary.
      */
      double rSpaceEnergy();

      /**
      * Get total Coulomb energy.
      *
      * Invokes kSpaceEnergy() and rSpaceEnergy() internally, and
      * thus recomputes k-space and r-space components as needed.
      */
      double energy();

      //@}
      /// \name Stress
      //@{

      /**
      * Compute kspace part of Coulomb stress.
      */
      virtual void computeStress() = 0;
   
      /**
      * Unset k-space stress.
      */
      virtual void unsetStress();

      /**
      * Get long-range k-space part of Coulomb stress.
      *
      * Recomputes iff necessary.
      */
      Tensor kSpaceStress();

      /**
      * Return short-range r-space part of Coulomb stress.
      * 
      * Directs associated pair potential to recompute if necessary.
      */
      Tensor rSpaceStress();

      /**
      * Get total Coulomb stress.  
      * 
      * Calls kSpaceStress() and rSpaceStress() internally, and thus
      * recomputes k-space and r-space components as needed.
      */
      Tensor stress();

      /**
      * Get total Coulomb stress.
      *
      * Equivalent to Tensor stress().
      */
      void computeStress(Tensor& stress);

      /**
      * Get diagonal components of Coulomb stress.
      */
      void computeStress(Vector& pressures);

      /**
      * Get total Coulomb pressure.
      */
      double pressure();

      /**
      * Get the total Coulomb pressure.
      * 
      * Equivalent to double pressure().
      */
      void computeStress(double& pressure);

      //@}

   protected:

      /// Short-range real space energy and stress contributions
      EwaldRSpaceAccumulator rSpaceAccumulator_;

      /// K-space part of Coulomb energy
      Setable<double> kSpaceEnergy_;

      /// K-space part of Coulomb stress.
      Setable<Tensor> kSpaceStress_;

      /// Have parameters been set?
      bool isInitialized_;

      /// Are waves and k-space potential up to date?
      /// Unset if boundary or parameters change
      bool hasWaves_;

   };

   /*
   * Are wavevectors and k-space potential up to date?
   */
   inline
   bool MdCoulombPotential::hasWaves()
   {  return hasWaves_; }

} 
#endif
