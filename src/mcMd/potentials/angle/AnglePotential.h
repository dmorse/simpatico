#ifndef MCMD_ANGLE_POTENTIAL_H
#define MCMD_ANGLE_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
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
   class AnglePotential : public ParamComposite
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
      * Default method throws an exception.This allows testing of subclasses 
      * that only work for MD simulation, and crash gracefully if used for MC.
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

      /**
      * Calculate the total angle potential energy for the System.
      */
      virtual double energy() const = 0;

      /**
      * Compute total angle potential isotropic pressure contribution.
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const = 0;

      /**
      * Compute x, y, z angle potential pressure components.
      *
      * \param stress (output) pressure components.
      */
      virtual void computeStress(Util::Vector& stress) const = 0;

      /**
      * Compute angle potential stress tensor contribution.
      *
      * \param stress (output) stress tensor contribution.
      */
      virtual void computeStress(Util::Tensor& stress) const = 0;
    
      //@}
   };

}

#endif
