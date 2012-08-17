#ifndef DDMD_POTENTIAL_CPP
#define DDMD_POTENTIAL_CPP

#include "Potential.h"

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   {}

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
   * Set a value for the total energy.
   */
   double Potential::energy()
   {  return energy_.value(); }

   /*
   * Set a value for the total energy.
   */
   void Potential::setEnergy(double energy)
   {  energy_ = energy; }

   /*
   * Mark the energy as unknown. 
   */
   void Potential::unsetEnergy()
   {  energy_.setNull(); }

   /*
   * Set a value for the total stress.
   */
   Tensor Potential::stress()
   {  return stress_.value(); }

   /*
   * Set a value for the total stress.
   */
   void Potential::setStress(const Tensor& stress)
   {  stress_ = stress; }

   /*
   * Mark the stress as unknown. 
   */
   void Potential::unsetStress()
   {  stress_.setNull(); }

}
#endif
