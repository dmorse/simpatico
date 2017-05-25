#include "Potential.h"

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   Potential::Potential()
    : stress_(),
      energy_(),
      reverseUpdateFlag_(false)
   { setClassName("Potential"); }

   /*
   * Destructor.
   */
   Potential::~Potential()
   {}

   /*
   * Set flag to identify if reverse communication is enabled.
   */
   void Potential::setReverseUpdateFlag(bool reverseUpdateFlag)
   { reverseUpdateFlag_ = reverseUpdateFlag; }

   /*
   * Get the value of the total energy.
   */
   double Potential::energy() const
   {  return energy_.value(); }

   /*
   * Set a value for the total energy.
   */
   void Potential::setEnergy(double energy)
   {  energy_.set(energy); }

   /*
   * Mark the energy as unknown. 
   */
   void Potential::unsetEnergy()
   {  energy_.unset(); }

   /*
   * Is the energy set?
   */
   bool Potential::isEnergySet() const
   {  return energy_.isSet(); }

   /*
   * Get the value for the total stress.
   */
   Tensor Potential::stress() const
   {  return stress_.value(); }

   /*
   * Set a value for the pressure.
   */
   double Potential::pressure() const
   {  
      double p = 0.0;
      for (int i = 0; i < Dimension; ++i) {
         p += stress_.value()(i, i);
      }
      p /= 3.0;
      return p;
   }

   /*
   * Set a value for the total stress.
   */
   void Potential::setStress(const Tensor& stress)
   {  stress_.set(stress); }

   /*
   * Mark the stress as unknown. 
   */
   void Potential::unsetStress()
   {  stress_.unset(); }

   /*
   * Is the stress set?
   */
   bool Potential::isStressSet() const
   {  return stress_.isSet(); }

   /*
   * Compute atomic forces and stress on all processors.
   * 
   * Default implementation just calls computeForces and computeStress.
   */
   #ifdef UTIL_MPI
   void Potential::computeForcesAndStress(MPI::Intracomm& communicator)
   #else
   void Potential::computeForcesAndStress();
   #endif
   { 
      computeForces(); 
      #ifdef UTIL_MPI
      computeStress(communicator); 
      #else
      computeForces();
      #endif
   }

   /*
   * Reduce energy from all processors.
   */
   #ifdef UTIL_MPI
   void Potential::reduceEnergy(double localEnergy, MPI::Intracomm& communicator)
   #else
   void Potential::reduceEnergy(double localEnergy)
   #endif
   {
      #ifdef UTIL_MPI
      double totalEnergy = 0.0; 
      communicator.Reduce(&localEnergy, &totalEnergy, 1, 
                          MPI::DOUBLE, MPI::SUM, 0);
      if (communicator.Get_rank() != 0) {
         totalEnergy = 0.0;
      }
      setEnergy(totalEnergy);
      #else
      setEnergy(localEnergy);
      #endif
   }

   /*
   * Reduce stress from all processors.
   */
   #ifdef UTIL_MPI
   void Potential::reduceStress(Tensor& localStress, MPI::Intracomm& communicator)
   #else
   void PairPotential::reduceStress(Tensor& localStress)
   #endif
   {
      #ifdef UTIL_MPI
      Tensor totalStress;
      communicator.Reduce(&localStress(0,0), &totalStress(0,0), 
                          Dimension*Dimension, MPI::DOUBLE, MPI::SUM, 0);
      if (communicator.Get_rank() != 0) {
         totalStress.zero();
      }
      setStress(totalStress);
      #else
      setStress(localStress);
      #endif
   }

   #ifdef UTIL_MPI
   /*
   * Is the potential in a valid internal state?
   */
   bool Potential::isValid(MPI::Intracomm& communicator) const
   {
      energy_.isValid(communicator);
      stress_.isValid(communicator);
      return true;
   }
   #endif

}
