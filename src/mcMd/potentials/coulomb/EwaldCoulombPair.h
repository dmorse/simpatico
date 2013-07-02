#ifndef MCMD_EWALD_COULOMB_PAIR_H
#define MCMD_EWALD_COULOMB_PAIR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>      // Member template
#include <mcMd/chemistry/AtomType.h>    // Member template parameter
#include <util/global.h>

#include <cmath>

namespace McMd
{

   using namespace Util;

   /**
   * Pair interaction for short-range part of Coulomb potential.
   */
   class EwaldCoulombPair
   {
   
   public:
   
      /**  
      * Set AtomTypes.
      */
      void setAtomTypes(Array<AtomType> atomTypes);

      /**
      * Returns interaction energy for a single pair of particles. 
      *
      * \param rsq square of distance between particles
      * \param i   type of particle 1
      * \param j   type of particle 2
      * \return    pair interaction energy
      */
      double energy(double rsq, int i, int j) const;
   
      /**
      * Returns ratio of scalar pair interaction force to pair separation.
      *
      * Multiply this quantity by the components of the separation vector
      * to obtain the force vector. A positive value for the return value
      * represents a repulsive force between a pair of particles.
      *
      * Precondition: The square separation rsq must be less than rCutoffSq.
      * If rsq > rCutoffSq, the return value is undefined (i.e., wrong).
      * Usage: Test for rsq < rCutoffSq before calling this function
      * \code
      * if (rsq < interaction.rCutoffSq(i, j)) {
      *    f = forceOverR(rsq, i, j);
      *    .....
      * }
      * \endcode
      *
      * \param rsq square of distance between particles
      * \param i type of particle 1
      * \param j type of particle 2
      * \return  force divided by distance 
      */
      double forceOverR(double rsq, int i, int j) const;
   
      //@}

   private:
   
      /// Maximum allowed value for nAtomType (# of particle types)
      static const int MaxAtomType = 2;
   
      /// Prefactor of qi*qj/r in Coulomb potential
      double c_;

      /// Prefactor 2*alpha/\sqrt{pi} in derivative of erfc
      double d_;

      /// Ewald inverse range parameter
      double alpha_;  
 
      /// Real space cutoff squared
      double rCutoffSq_;  

      /// Pointer to array of AtomTypes
      const Array<AtomType>* atomTypesPtr_;

      /**
      * Set all private member constants. Called by CoulombPotential.
      */
      void set(double epsilon, double alpha, double rCutoff);

   //friends:

      friend class CoulombPotential;

   };
  
   // inline methods 
 
   /* 
   * Calculate interaction energy for a pair, as function of squared distance.
   */
   inline 
   double EwaldCoulombPair::energy(double rsq, int i, int j) const 
   {
      if (rsq < rCutoffSq_) {
         double qi = (*atomTypesPtr_)[i].charge();
         double qj = (*atomTypesPtr_)[j].charge();
         double r = sqrt(rsq);
         return c_*qi*qj*erfc(alpha_*r)/r;
      } else {
         return 0.0;
      }
   }
   
   /* 
   * Calculate force/distance for a pair as function of squared distance.
   */
   inline 
   double EwaldCoulombPair::forceOverR(double rsq, int i, int j) const
   {
      double r = sqrt(rsq);
      double qi = (*atomTypesPtr_)[i].charge();
      double qj = (*atomTypesPtr_)[j].charge();
      return c_*qi*qj*( erfc(alpha_*r) + d_*exp(-alpha_*rsq) )/(r*rsq);
   }

}
#endif
