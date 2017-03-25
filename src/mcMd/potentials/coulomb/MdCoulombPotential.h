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
   * Coulomb potential base class.
   *
   * This class manages parameters for an Ewald summation 
   * algorithm for calculating Coulomb energies and forces.
   * Functions that compute different contributions to the
   * energies and forces are divided among several classes.
   * The short-range pair potential and pair force for a 
   * single pair can be calculated by member functions of 
   * an associated instance of the EwaldCoulombPair class.
   * This class implements the calculation of the k-space 
   * (Fourier space) contribution to the total energy. 
   *
   * An association between this class and an EwaldCoulombPair
   * is created by the setPairInteraction method, which sets
   * a private pointer to the EwaldCoulombPair in this class. 
   * Functions in this class that set or modify the epsilon, 
   * alpha, and rCutoff member variables also automatically
   * set or modify all related constants in the associated
   * EwaldCoulombPair. By design, this is the only possible 
   * way to set or modify the parameters of the EwaldCoulombPair, 
   * and is guarantees that the parameters used by these two
   * classes are always consistent.
   *
   * Implementation notes: This class is a friend of the
   * EwaldCoulombPair class, and calls the private member 
   * function EwaldCoulombPair::set to set parameters of 
   * the associated EwaldCoulombPair. The EwaldCoulombPair 
   * class does not provide any public member functions 
   * that can set these parameters.
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

      /// \name Initialization
      //@{

      /**
      * Read parameters and initialize.
      *
      * \param in input stream
      */
      //virtual void readParameters(std::istream& in);

      /**
      * Generate wavevectors for the current boundary.
      */
      virtual void makeWaves() = 0;

      //@}
      /// \name System energy and stress.
      //@{

      /**
      * Add k-space Coulomb forces for all atoms.
      */
      virtual void addForces() = 0;

      /**
      * Calculate the long range kspace part of Coulomb energy.
      */
      virtual void computeEnergy() = 0;

      // old convention for test. huada
      //virtual void computeKSpaceStress(double& stress) = 0;
      //virtual void computeKSpaceStress(Util::Vector& stress) = 0;
      //virtual void computeKSpaceStress(Util::Tensor& stress) = 0;

      /**
      * Compute kspace part of Coulomb stress and pressure at the same time.
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
      //virtual double kSpacePressure() const = 0;
      virtual Tensor kSpaceStress() const = 0;

      // Return R-space contributions
      virtual double rSpaceEnergy() const = 0;
      //virtual double rSpacePressure() const = 0;
      virtual Tensor rSpaceStress() const = 0;

      // Return total energy and stress
      virtual double energy() = 0;
      //virtual double pressure() const = 0;
      virtual Tensor stress() = 0;

      //@}

   protected:

      /// Has readParam been called?
      bool isInitialized_;

   };
} 
#endif
