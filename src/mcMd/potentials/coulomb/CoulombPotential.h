#ifndef MCMD_COULOMB_POTENTIAL_H
#define MCMD_COULOMB_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <mcMd/chemistry/AtomType.h>     // Member template parameter
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
   class CoulombPotential : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      CoulombPotential(System& system);

      /**
      * Destructor (does nothing).
      */
      virtual ~CoulombPotential();

      /// \name Initialization
      //@{

      /**
      * Create association with an EwaldCoulombPair pair interaction.
      *
      * \param pairInteraction associated EwaldCoulombPair object
      */
      void setPairInteraction(EwaldCoulombPair& pairInteraction);

      /**
      * Read parameters and initialize.
      *
      * \param in input stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Generate wavevectors for the current boundary.
      */
      virtual void makeWaves() = 0;

      //@}
      /// \name System energy and stress.
      //@{

      /**
      * Calculate the long range kspace part of Coulomb energy.
      */
      virtual double kspaceEnergy() = 0;

      /**
      * Add k-space Coulomb forces for all atoms.
      */
      virtual void addKSpaceForces() = 0;

      /**
      * Compute kspace part of Coulomb pressure.
      *
      * \param stress (output) pressure
      */
      virtual void computeKSpaceStress(double& stress) = 0;

      /**
      * Compute kspace part of Coulomb x, y, z pressures.
      *
      * \param stress (output) pressures
      */
      virtual void computeKSpaceStress(Util::Vector& stress) = 0;

      /**
      * Compute kspace part of Coulomb stress tensor.
      *
      * \param stress (output) stress tensor
      */
      virtual void computeKSpaceStress(Util::Tensor& stress) = 0;
    
      //@}
      /// \name Accessors (const)
      //@{

      /**
      * Get value of dielectric constant.
      */
      double epsilon() const;

      /**
      * Get value of Ewald inverse range parameter.
      */
      double alpha() const;

      /**
      * Get cutoff distance for short range interaction.
      */
      double rCutoff() const;

      /**
      * Current number of wavevectors.
      */
      virtual int nWave() const = 0;

      //@}

   protected:

      /// Pointer to associated Simulation
      Simulation* simulationPtr_;

      /// Pointer to associated System
      System* systemPtr_;

      /// Pointer to associated Boundary
      Boundary* boundaryPtr_;

      /// Pointer to array of AtomTypes.
      const Array<AtomType>* atomTypesPtr_;

      /// Pointer to associated short range pair interaction.
      EwaldCoulombPair* pairInteractionPtr_;

      /// Dielectric constant / permittivity.
      double epsilon_;

      /// Inverse range parameter.
      double alpha_; 

      /// Cutoff distance for short range pair interaction.
      double rCutoff_;

      /// Has readParam been called?
      bool isInitialized_;

   };

   // Inline functions

   /*
   * Get value of epsilon (dielectric) parameter.
   */
   inline double CoulombPotential::epsilon() const
   {  return epsilon_; }

   /*
   * Get value of inverse range parameter.
   */
   inline double CoulombPotential::alpha() const
   {  return alpha_; }

   /*
   * Get cutoff distance for short range interaction.
   */
   inline double CoulombPotential::rCutoff() const
   {  return rCutoff_; }

} 
#endif
