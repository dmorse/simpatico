#ifdef  MCMD_PERTURB
#ifndef MCMD_MC_ENERGY_PERTURBATION_H
#define MCMD_MC_ENERGY_PERTURBATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/perturb/LinearPerturbation.h>

namespace McMd
{

   class McSystem;

   /**
   * A Perturbation with a variable inverse temperature.
   *
   * An McEnergyPerturbation is a series of systems with the same Hamiltonian
   * at different temperatures. The perturbation parameter is beta = 1/kT.
   *
   * \ingroup McMd_Perturb_Module
   */
   class McEnergyPerturbation : public LinearPerturbation<McSystem>
   {

   public:
  
      /**
      * Constructor. 
      */
      McEnergyPerturbation(McSystem& system, int size, int rank);
  
      /**
      * Read beta parameter (inverse temperature) from file.
      *
      * \param in input stream (file or std in).
      */ 
      virtual void readParameters(std::istream& in);
 
      /**
      * Set inverse temperature of the parent system.
      */
      virtual void setParameter();

      /**
      * Get inverse temperature of the parent system.
      *
      * \param i index of the perturbation parameter. 
      * In this case, i is always 0 since nParameters = 1.
      */
      virtual double parameter(int i) const;

      /**
      * Get derivative of the Boltzmann weight with respect to the perturbation 
      * parameter.
      *
      * For a system in which W = beta*H depends on a parameter x, this method
      * returns the value of dW/dx for the current system configuration.
      *
      * \param i index of the perturbation parameter. 
      * In this case, i is always 0 since nParameters = 1.
      */
      virtual double derivative(int i) const;

   private:  

      /*
      Number of perturbation parameters associated with a System.
      nParameters = 1 for McEnergyPerturbation.
      */
      int nParameters_;

   };

}
#endif  // ifndef ENERGY_PERTURBATION_H
#endif  // ifdef MCMD_PERTURB
