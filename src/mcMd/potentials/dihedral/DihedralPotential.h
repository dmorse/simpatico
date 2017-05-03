#ifndef MCMD_DIHEDRAL_POTENTIAL_H
#define MCMD_DIHEDRAL_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <mcMd/potentials/misc/EnergyCalculator.h>   // base class
#include <mcMd/potentials/misc/StressCalculator.h>   // base class

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
   * Interface for a Dihedral Potential.
   *
   * \ingroup McMd_Dihedral_Module
   */
   class DihedralPotential : public ParamComposite,
                             public EnergyCalculator, public StressCalculator
   
   {

   public:

      /**
      * Constructor.
      */
      DihedralPotential();

      /**
      * Destructor (does nothing)
      */
      virtual ~DihedralPotential();

      /// \name Dihedral interaction interface
      //@{ 

      /**
      * Returns potential energy for one dihedral.
      *
      *     0   2    3
      *     o   o----o
      *      \ /
      *       o 
      *       1 
      *
      * \param R1     bond vector r1 - r0 from atom 0 to 1
      * \param R2     bond vector r2 - r1 from atom 1 to 2
      * \param R3     bond vector r3 - r2 from atom 2 to 3
      * \param type   type index for dihedral
      */
      virtual
      double energy(const Vector& R1, const Vector& R2, const Vector& R3,
                    int type) const = 0;
 
      // Prevent hiding of inherited function energy();
      using EnergyCalculator::energy;
    
      /**
      * Returns derivatives of energy with respect to bond vectors.
      *
      * \param R1    bond vector from atom 0 to 1 (input)
      * \param R2    bond vector from atom 1 to 2 (input)
      * \param R3    bond vector from atom 2 to 3 (input)
      * \param F1    derivative of energy w/respect to R1 (output)
      * \param F2    derivative of energy w/respect to R2 (output)
      * \param F3    derivative of energy w/respect to R3 (output)
      * \param type  type index for dihedral
      */
      virtual
      void force(const Vector& R1, const Vector& R2, const Vector& R3,
                 Vector& F1, Vector& F2, Vector& F3, int type) const = 0;

      /**
      * Modify an interaction parameter, identified by a string.
      *
      * \param name  parameter name
      * \param type  dihedral type index 
      * \param value new value of parameter
      */
      virtual void set(std::string name, int type, double value) = 0;

      /**
      * Get an interaction parameter value, identified by a string.
      *
      * \param  name  parameter name
      * \param  type  dihedral type index
      * \return value of parameter
      */
      virtual double get(std::string name, int type) const = 0;

      /**
      * Return name of pair interaction class (e.g., "HarmonicDihedral").
      */
      virtual std::string interactionClassName() const = 0;

      //@}
      /// \name System energy, force, and stress
      //@{ 

      /**
      * Compute and return the dihedral potential energy for one Atom.
      *
      * Default method throws an exception.This allows testing of subclasses 
      * that only work for MD simulation, and crash gracefully if used for MC.

      * \param  atom Atom object of interest
      * \return bond potential energy of atom
      */
      virtual double atomEnergy(const Atom& atom) const
      {  
         UTIL_THROW("Unimplemented method"); 
         return 0.0; // Never reached, but avoids compiler warning.
      }

      /**
      * Add dihedral potential forces to all atomic forces.
      *
      * Default version throws an exception.This allows testing of subclasses 
      * that only work for MC simulation, and crash gracefully if used for MD.
      */
      virtual void addForces()
      {  UTIL_THROW("Unimplemented method"); }
    
   };

} 
#endif
