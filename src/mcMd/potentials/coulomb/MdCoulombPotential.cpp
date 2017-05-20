#ifndef MD_COULOMB_POTENTIAL_CPP
#define MD_COULOMB_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdCoulombPotential.h" 
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MdCoulombPotential::MdCoulombPotential()
    : isInitialized_(false),
      hasWaves_(false)
   {  setClassName("CoulombPotential"); }

   /*
   * Destructor (does nothing)
   */
   MdCoulombPotential::~MdCoulombPotential()
   {}

   /*
   * Modify an interaction parameter, identified by a string.
   *
   * This default implementation throw an Exception, but is 
   * called only if virtual function is not re-implemented by 
   * derived class.
   */
   void MdCoulombPotential::set(std::string name, double value)
   {
      UTIL_THROW("Unimplemented virtual set method");  
   };

   /*
   * Get an interaction parameter value, identified by a string.
   *
   * This default implementation throw an Exception, but is 
   * called only if virtual function is not re-implemented by 
   * derived class.
   */
   double MdCoulombPotential::get(std::string name) const 
   {
      UTIL_THROW("Unimplemented virtual get method");  
      return 0.0; 
   };

   /*
   * Unset precomputed wavevectors, influence function, energy and stress.
   */
   void MdCoulombPotential::unsetWaves()
   {  
      hasWaves_ = false; 
      unsetEnergy();
      unsetStress();
   }

   /*
   * Unset k-space energy (does not unset rSpace or waves).
   */
   void MdCoulombPotential::unsetEnergy()
   {  kSpaceEnergy_.unset(); }

   /*
   * Get kspace energy (compute iff necessary).
   */
   double MdCoulombPotential::kSpaceEnergy()
   {
      if (!kSpaceEnergy_.isSet()) {
         computeEnergy();
      }
      return kSpaceEnergy_.value();
   }

   /*
   * Get rspace energy (compute iff necessary).
   */
   double MdCoulombPotential::rSpaceEnergy()
   {  return rSpaceAccumulator_.rSpaceEnergy(); }

   /*
   * Get total Coulomb energy (recompute as needed). 
   */
   double MdCoulombPotential::energy()
   {
      double temp;
      temp  = kSpaceEnergy();
      temp += rSpaceEnergy();
      return temp;
   }

   /*
   * Unset k-space stress contribution (does not unset rSpace).
   */
   void MdCoulombPotential::unsetStress()
   {  kSpaceStress_.unset(); }

   /*
   * Get k-space stress contribution.
   */
   Tensor MdCoulombPotential::kSpaceStress() 
   {
      if (!kSpaceStress_.isSet()) {
         computeStress();
      }
      return kSpaceStress_.value();
   }

   /*
   * Get r-space stress contribution (recompute as needed).
   */
   Tensor MdCoulombPotential::rSpaceStress()
   {  return rSpaceAccumulator_.rSpaceStress(); }

   /*
   * Get total Coulomb stress (recompute as needed).
   */
   Tensor MdCoulombPotential::stress()
   {
      Tensor temp;
      temp = kSpaceStress();
      temp += rSpaceAccumulator_.rSpaceStress();
      return temp;
   }
 
   /*
   * Get total Coulomb stress (recompute as needed).
   */
   double MdCoulombPotential::pressure()
   {
      Tensor temp = stress();
      double value = 0.0;
      for (int i = 0; i < Dimension; ++i) {
         value += temp(i, i); 
      }
      return value/3.0;
   }
 
   /*
   * Get total Coulomb stress.
   */
   void MdCoulombPotential::computeStress(Tensor& stress)
   {  
      stress = kSpaceStress();
      stress+= rSpaceAccumulator_.rSpaceStress();
   }

   /**
   * Get diagonal components of Coulomb stress.
   */
   void MdCoulombPotential::computeStress(Vector& pressures)
   {
      Tensor temp = stress();
      for (int i = 0; i < Dimension; ++i) {
         pressures[i] = temp(i,i);
      }
   }

   /**
   * Get Coulomb pressure.
   * 
   * Equivalent to double pressure().
   */
   void MdCoulombPotential::computeStress(double& pressure)
   {
      Tensor temp = stress();
      pressure = 0.0;
      for (int i = 0; i < Dimension; ++i) {
         pressure += temp(i, i); 
      }
      pressure /= 3.0;
   }

} 
#endif
