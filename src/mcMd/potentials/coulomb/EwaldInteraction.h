#ifndef MCMD_EWALD_INTERACTION_H
#define MCMD_EWALD_INTERACTION_H

/*
 * Simpatico - Simulation Package for Polymeric and Molecular Liquids
 *
 * Copyright 2010 - 2012, David Morse (morse012@umn.edu)
 * Distributed under the terms of the GNU General Public License.
 */

#include <util/param/ParamComposite.h>
#include <util/global.h>

#include <math.h>

namespace McMd
{

   using namespace Util;

   /**
   * Implementation of pairwise r- and k-space Ewald interaction.
   *
   * This class defines the standard Ewald potential.
   *
   * \ingroup McMd_Coulomb_Module
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
      * \param other EwaldInteraction to be assigned.
      */
      EwaldInteraction& operator = (const EwaldInteraction& other);


      /**
      * Default destructor.
      */
      ~EwaldInteraction() {}

      /// \name Mutators
      //@{ 

      /**
      * Read epsilon, alpha, rCutoff, and kCutoff.
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
      /// \name Accessors
      //@{ 

      /**
      * Returns interaction energy for a single pair of atoms. 
      *
      * \param rsq  square of distance between atoms
      * \return  pair interaction energy
      */
      double rSpaceEnergy(double rSq, double qProduct) const;

      /**
      * Returns interaction energy for a single pair of atoms. 
      *
      * \param rsq square of distance between atoms
      * \return    pair interaction energy
      */
      double kSpacePotential(double rhoSq, double prefactor) const;

      /**
      * Returns ratio of scalar pair interaction force to pair separation.
      *
      * Multiply this quantity by the components of the separation vector
      * to obtain the force vector. A positive value for the return value
      * represents a repulsive force between a pair of atoms.
      *
      * Precondition: The distance squared rsq must be less than cutoffSq.
      * If rsq > cutoffSq, the return value is undefined (i.e., invalid).
      * Usage: Test for rsq < cutoffSq before calling this function
      * \code
      * if (rsq < interaction.cutoffSq(i, j)) {
      *    f = forceOverR(rsq, i, j);
      *    .....
      * }
      * \endcode
      *
      * \param rsq square of distance between atoms
      * \return    force divided by distance 
      */
      double rSpaceForceOverR(double rSq, double qProduct) const;
  
      /**
      * Get real space cutoff squared.
      */
      double rSpaceCutoffSq() const;

      /**
      * Get real space cutoff.
      */
      double rSpaceCutoff() const;

      /**
      * Get Ewald paramter alpha.
      *
      * \return    alpha_
      */
      double alpha() const;
 
      /**
      * Get dielectric permittivity.
      *
      * \return    epsilon_
      */
      double epsilon() const;
 
      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      */
      double get(std::string name) const;

      //@}

   protected:

      /// Physical Parameters.
      double epsilon_;          ///< Dielectric permittivity.

      /// Algorithmic parameters for Ewald potential.
      double alpha_;            ///< alpha = (1 / (sigma*sqrt(2)) ).
      double rSpaceCutoff_;     ///< Ewald potential real space cutoff.
      double rSpaceCutoffSq_;   ///< Real space cutoff squared.

      /// prefactors for real space potential
      double ce_;
      double cf_;

      /**
      * Was this object initialized by calling (read|load)Parameters ?
      */
      bool  isInitialized;

   };

   // Inline methods 

   /* 
   * Return medium dielectric permittivity.
   */
   inline double EwaldInteraction::epsilon() const
   { return epsilon_; }

   /* 
   * Return Ewald mearing parameter alpha.
   */
   inline double EwaldInteraction::alpha() const
   { return alpha_; }

   /* 
   * Return real space cutoff distance in Ewald method.
   */
   inline double EwaldInteraction::rSpaceCutoff() const
   {  return rSpaceCutoff_; }

  /* 
   * Return real space cutoff distance squared in Ewald method.
   */
   inline double EwaldInteraction::rSpaceCutoffSq() const
   { return rSpaceCutoffSq_; }

   /* 
   * Calculate r-space energy for a pair of charges.
   */
   inline double EwaldInteraction::rSpaceEnergy(double rsq, double qProduct) const 
   {
      double r = sqrt(rsq);
      return ce_*qProduct*erfc(alpha_*r)/r;
   }

   /*
   * Calculate r-space force/distance for a pair of charges.
   */
   inline double EwaldInteraction::rSpaceForceOverR(double rSq, double qProduct) const 
   {
      double r = sqrt(rSq);
      double x = alpha_*r;
      return ce_*qProduct*(erfc(x) + cf_*r*exp(-x*x))/(r*rSq); 
   }

   /* 
   * Calculate k-space potential from wavenumber kSq.
   */
   inline double EwaldInteraction::kSpacePotential(double rhoSq, double prefactor) const
   {
      return rhoSq*prefactor;
   }

}
#endif
