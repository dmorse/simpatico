#ifndef MCMD_DIHEDRAL_POTENTIAL_H
#define MCMD_DIHEDRAL_POTENTIAL_H

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
   * Interface for a Dihedral Potential.
   *
   * \ingroup McMd_Dihedral_Module
   */
   class DihedralPotential : public ParamComposite
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
      *     1   3    4
      *     o   o----o
      *      \ /
      *       o 
      *       2 
      *
      * \param R1     bond vector from atom 1 to 2.
      * \param R2     bond vector from atom 2 to 3.
      * \param R3     bond vector from atom 3 to 4.
      * \param type   type of dihedral.
      */
      virtual
      double energy(const Vector& R1, const Vector& R2, const Vector& R3,
                    int type) const = 0;
 
      /**
      * Returns derivatives of energy with respect to bond vectors forming the
      * dihedral group.
      *
      * \param R1     bond vector from atom 1 to 2.
      * \param R2     bond vector from atom 2 to 3.
      * \param R3     bond vector from atom 3 to 4.
      * \param F1     return force along R1 direction.
      * \param F2     return force along R2 direction.
      * \param F3     return force along R2 direction.
      * \param type   type of dihedral.
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
      * \param  atom Atom object of interest
      * \return bond potential energy of atom
      */
      virtual double atomEnergy(const Atom& atom) const
      {  UTIL_THROW("Unimplemented method"); }

      /**
      * Compute and return the total Dihedral energy for the associated System.
      */
      virtual double energy() const = 0;

      /**
      * Add dihedral potential forces to all atomic forces.
      */
      virtual void addForces()
      {  UTIL_THROW("Unimplemented method"); }

      /**
      * Compute total dihedral pressure
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const = 0;

      /**
      * Compute x, y, z dihedral pressures.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const = 0;

      /**
      * Compute dihedral stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const = 0;
    
   };

} 
#endif
