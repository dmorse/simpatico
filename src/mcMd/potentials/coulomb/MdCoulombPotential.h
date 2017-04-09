#ifndef MD_COULOMB_POTENTIAL_H
#define MD_COULOMB_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class

#include <complex>

namespace Util
{
   class Vector;
   class Tensor;
}

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

      /// \name System energy and stress.
      //@{

      /**
      * Generate wavevectors for the current boundary.
      */
      virtual void makeWaves() = 0;

      /**
      * Add k-space Coulomb forces for all atoms.
      */
      virtual void addForces() = 0;

      /**
      * Unset the long range kspace part of Coulomb energy.
      */
      virtual void unsetEnergy() = 0;

      /**
      * Calculate the long range kspace part of Coulomb energy.
      */
      virtual void computeEnergy() = 0;

      /**
      * Unset the long range kspace part of Coulomb stress.
      */
      virtual void unsetStress() = 0;

      /**
      * Compute kspace part of Coulomb stress.
      *
      * \param stress (output) pressure
      */
      virtual void computeStress() = 0;
   
      //@}
      /// \name Accessors (const)
      //@{

      /**
      * Current number of wavevectors.
      */
      virtual int nWave() const = 0;

      // Return K-space contributions
      virtual double kSpaceEnergy() const = 0;
      virtual Tensor kSpaceStress() const = 0;

      // Return R-space contributions
      virtual double rSpaceEnergy() const = 0;
      virtual Tensor rSpaceStress() const = 0;

      // Return total energy and stress
      virtual double energy() = 0;
      virtual Tensor stress() = 0;

      //@}

   protected:

      /// Has readParam been called?
      bool isInitialized_;

   };
} 
#endif
