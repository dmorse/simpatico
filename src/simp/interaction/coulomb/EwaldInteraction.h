#ifndef SIMP_EWALD_INTERACTION_H
#define SIMP_EWALD_INTERACTION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/global.h>

#include <math.h>

namespace Simp
{

   using namespace Util;

   /**
   * Implementation of r-space and k-space Ewald Coulomb interactions.
   *
   * This class defines the standard Ewald decomposition of a coulomb
   * potential. The total pair potential is 1/(4 pi epsilon r). The
   * short-range part is erfc(alpha r)/(4 pi epsilon r). The Gaussian
   * smeared k-space potential for square wavenumber kSq is given by
   * V(k) = exp(-kSq/(4*alpha^2))/(epsilon kSq).
   *
   * \ingroup Simp_Coulomb_Module
   */
   class EwaldInteraction : public ParamComposite
   {
   public:

      /**
      * Default constructor.
      */
      EwaldInteraction();

      /**
      * Copy constructor.
      *
      * \param other EwaldInteraction to be copied.
      */
      EwaldInteraction(const EwaldInteraction& other);

      /**
      * Assignment.
      *
      * \param other EwaldInteraction from which data is copied.
      */
      EwaldInteraction& operator = (const EwaldInteraction& other);

      /**
      * Default destructor.
      */
      ~EwaldInteraction() 
      {}

      /// \name Mutators
      //@{ 

      /**
      * Read epsilon, alpha, and rCutoff.
      *
      * \param in  input parameter stream 
      */
      void readParameters(std::istream &in);

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

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param value  new value of parameter
      */
      void set(std::string name, double value);

      //@}
      /// \name Energy and Force Computations
      //@{

      /**
      * Returns r-space interaction energy for a single pair of atoms. 
      *
      * rSpaceEnergy = qProduct*erfc(r*alpha)/(4*pi*epsilon*r)
      *
      * \param rSq  square of distance between atoms
      * \param qProduct  product of charges
      * \return  short range part of pair energy 
      */
      double rSpaceEnergy(double rSq, double qProduct) const;

      /**
      * Return ratio of scalar pair interaction force to pair separation.
      *
      * Multiply this quantity by the components of the separation vector
      * to obtain the force vector. A positive value for the return value
      * represents a repulsive force between a pair of atoms.
      *
      * Precondition: The distance squared rSq must be less than cutoffSq.
      * If rSq > cutoffSq, the return value is undefined (i.e., invalid).
      * Usage: Test for rSq < cutoffSq before calling this function
      * \code
      * if (rSq < interaction.rSpaceCutoffSq()) {
      *    forceOverR = rSpaceForceOverR(rsq, i, j);
      *    .....
      * }
      * \endcode
      *
      * \param rSq square of distance between atoms
      * \param qProduct product of charges
      * \return force magnitude divided by distance 
      */
      double rSpaceForceOverR(double rSq, double qProduct) const;
  
      /**
      * Return regularized Fourier-space potential.
      *
      * \param kSq square of wavenumber
      * \return exp(-kSq/(4*alpha*alpha))/(epsilon*ksq)
      */
      double kSpacePotential(double kSq) const;

      //@}
      /// \name Parameter Accessors
      //@{ 

      /**
      * Get dielectric permittivity.
      *
      * \return epsilon
      */
      double epsilon() const;
 
      /**
      * Get Ewald smearing parameter alpha (inverse length).
      *
      * \return alpha
      */
      double alpha() const;
 
      /**
      * Get real space cutoff.
      */
      double rSpaceCutoff() const;

      /**
      * Get real space cutoff squared.
      */
      double rSpaceCutoffSq() const;

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name parameter name
      */
      double get(std::string name) const;

      //@}

   private:

      /// Dielectric permittivity (arbitrary units).
      double epsilon_;          

      /// Ewald smearing parameter. Inverse length: 1/alpha = sigma*sqrt(2).
      double alpha_; 

      /// Real space cutoff for short range pair potential.
      double rSpaceCutoff_;  

      // Derived constants
      double rSpaceCutoffSq_; 
      double ce_;
      double cf_;
      double cg_;

      /**
      * Was this object initialized by calling (read|load)Parameters ?
      */
      bool isInitialized_;

      /// Compute and set values of all derived constants
      void setDerivedConstants();

   };

   // Inline methods 

   /* 
   * Return medium dielectric permittivity.
   */
   inline double EwaldInteraction::epsilon() const
   {  return epsilon_; }

   /* 
   * Return Ewald smearing parameter alpha.
   */
   inline double EwaldInteraction::alpha() const
   {  return alpha_; }

   /* 
   * Return real space cutoff distance.
   */
   inline 
   double EwaldInteraction::rSpaceCutoff() const
   {  return rSpaceCutoff_; }

  /* 
   * Return real space cutoff distance squared.
   */
   inline 
   double EwaldInteraction::rSpaceCutoffSq() const
   {  return rSpaceCutoffSq_; }

   /* 
   * Compute and return r-space energy for a pair of charges.
   */
   inline 
   double EwaldInteraction::rSpaceEnergy(double rSq, double qProduct) 
   const 
   {
      double r = sqrt(rSq);
      return ce_*qProduct*erfc(alpha_*r)/r;
   }

   /*
   * Compute and return (r-space force) / distance for a pair of charges.
   */
   inline 
   double EwaldInteraction::rSpaceForceOverR(double rSq, double qProduct) 
   const 
   {
      double r = sqrt(rSq);
      double x = alpha_*r;
      return ce_*qProduct*(erfc(x) + cf_*r*exp(-x*x))/(r*rSq); 
   }

   /* 
   * Calculate k-space potential from squared wavenumber kSq.
   */
   inline
   double EwaldInteraction::kSpacePotential(double kSq) const
   {  return exp(cg_*kSq)/(kSq*epsilon_); }

}
#endif
