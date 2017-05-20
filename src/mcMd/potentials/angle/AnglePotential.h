#ifndef MCMD_ANGLE_POTENTIAL_H
#define MCMD_ANGLE_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <mcMd/potentials/misc/EnergyCalculator.h>   // base class
#include <mcMd/potentials/misc/StressCalculator.h>   // base class

#include <util/random/Random.h>
#include <string>

namespace Util
{
   class Vector;
   class Tensor;
}

namespace McMd
{

   using namespace Util;

   class Atom;

   /**
   * Interface for a Angle Interaction.
   *
   * \ingroup McMd_Angle_Module
   */
   class AnglePotential : public ParamComposite,
                         public EnergyCalculator, public StressCalculator
   
   {

   public:

      /**
      * Constructor.
      */
      AnglePotential();

      /**
      * Destructor (does nothing)
      */
      virtual ~AnglePotential();

      /// \name Angle Interaction interface
      //@{

      /**
      * Returns potential energy for one angle.
      *
      * \param cosTheta  cosine of the bend angle.
      * \param type      type of bend angle.
      */
      virtual double energy(double cosTheta, int type) const = 0;
 
      // Prevent hiding of inherited function energy();
      using EnergyCalculator::energy;
    
      /**
      * Returns forces along two bonds at the angle, for use in MD and stress
      * calculation.
      *
      * \param R1     bond vector from atom 1 to 2.
      * \param R2     bond vector from atom 2 to 3.
      * \param F1     return force along R1 direction.
      * \param F2     return force along R2 direction.
      * \param type   type of angle.
      */
      virtual void force(const Vector& R1, const Vector& R2,
                         Vector& F1, Vector& F2, int type) const = 0;

      /**
      * Return bond angle chosen from equilibrium distribution.
      *
      * This function returns a bond length chosen from the Boltzmann
      * distribution of lengths for bonds of random orientation. The
      * distribution P(l) of values of the length l is proportional 
      * to l*l*exp[-beta*phi(l) ], where phi(l) is the bond energy. 
      *
      * \param random  pointer to random number generator object.
      * \param beta  inverse temperature
      * \param type  integer angle typd index
      * \return  angle chosen from equilibrium distribution.
      */
      virtual 
      double randomAngle(Util::Random* random, double beta, int type) 
             const = 0;

      /**
      * Return bond angle chosen from equilibrium distribution.
      *
      * This function returns a bond length chosen from the Boltzmann
      * distribution of lengths for bonds of random orientation. The
      * distribution P(l) of values of the length l is proportional 
      * to l*l*exp[-beta*phi(l) ], where phi(l) is the bond energy. 
      *
      * \param random  pointer to random number generator object.
      * \param beta  inverse temperature
      * \param type  bond type
      * \return  angle cosine chosen from equilibrium distribution.
      */
      virtual 
      double randomCosineAngle(Util::Random* random, double beta, int type) 
             const = 0;

      /**
      * Modify an interaction parameter, identified by a string.
      *
      * \param name  parameter name
      * \param type  angle type index 
      * \param value new value of parameter
      */
      virtual void set(std::string name, int type, double value) = 0;

      /**
      * Get an interaction parameter value, identified by a string.
      *
      * \param name parameter name
      * \param type angle type index 1
      */
      virtual double get(std::string name, int type) const = 0;

      /**
      * Return name of pair interaction class (e.g., "HarmonicAngle").
      */
      virtual std::string interactionClassName() const = 0;

      //@}
      /// \name System energy, force, and stress.
      //@{

      /**
      * Calculate the covalent bond energy for one Atom.
      *
      * Default implementation throws an exception. This allows testing of 
      * subclasses that only work for MD simulation, and crash gracefully 
      * if used for MC.
      *
      * \param  atom Atom object of interest
      * \return bond potential energy of atom
      */
      virtual double atomEnergy(const Atom& atom) const
      {  
         UTIL_THROW("Unimplemented method"); 
         return 0.0; // Never reached, but avoids compiler warning.
      }

      /**
      * Add bond forces to all atomic forces.
      *
      * Default version throws an exception.This allows testing of subclasses 
      * that only work for MC simulation, and crash gracefully if used for MD.
      */
      virtual void addForces()
      {  UTIL_THROW("Unimplemented method"); }

      //@}
   };

}
#endif
