#ifdef  MCMD_PERTURB
#ifndef MCMD_LINEAR_PERTURBATION_H
#define MCMD_LINEAR_PERTURBATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Perturbation.h"

namespace McMd
{

   /**
   * A Perturbation that is a linear function of a parameter.
   * 
   * A LinearPerturbation is a Perturbation for which the associated 
   * Boltzmann weight W(X, p) = H/kT for a system with a microstate X 
   * and a perturbation parameter p is a linear function of a 
   * parameter x. 
   *
   * One example of this is a perturbation in which the parameter 
   * p is equal to the inverse temperature beta = 1/kT. Another is 
   * when the Hamiltonian depends linearly upon some parameter p.
   */
   template <class SystemType>
   class LinearPerturbation : public Perturbation
   {
   
   public:

      using Perturbation::getNParameters;
      using Perturbation::parameter;
      using Perturbation::derivative;

      /**
      * Constructor. 
      *
      * Retains a pointer to the parent system, and invokes the
      * system.setPerturbation(*this) method to associate the new 
      * object with the SystemType object system.
      *
      * \param system parent system.
      * \param size   number of systems (communicator size)
      * \param rank   id of this system (communicator rank)
      */
      LinearPerturbation(SystemType& system, int size, int rank);

      /**
      * Destructor.
      */
      virtual ~LinearPerturbation();

      /**
      * Returns the difference W(partner) - W(current).
      */
      virtual double difference(DArray<double> iPartnerParameter) const;

      /**
      * Get the associated System by reference.
      */
      SystemType& system() const;

   private:

      /// Pointer to associated System.
      SystemType* systemPtr_;

   };

   // Methods

   /*
   * Constructor. 
   */
   template <class SystemType>
   LinearPerturbation<SystemType>
        ::LinearPerturbation(SystemType& system, int size, int rank)
    : Perturbation(size, rank),
      systemPtr_(&system)
   {}

   template <class SystemType>
   LinearPerturbation<SystemType>::~LinearPerturbation()
   {}

   /*
   * Returns the difference W(partner) - W(current).
   */
   template <class SystemType>
   double LinearPerturbation<SystemType>::difference(DArray<double> iPartnerParameter) const
   {  
      int nParameters = getNParameters();
      double difference = 0.0;
      for (int i = 0; i < nParameters; ++i) {
         difference += (iPartnerParameter[i] - parameter(i)) * derivative(i); 
      }
      return difference;
   }

   /*
   * Return associated System by reference.
   */
   template <class SystemType>
   inline SystemType& LinearPerturbation<SystemType>::system() const
   {  return *systemPtr_; }

}
      
#endif  // ifndef LINEAR_PERTURBATION_H
#endif  // ifdef  MCMD_PERTURB
