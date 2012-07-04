#ifndef DDMD_ANGLE_POTENTIAL_H
#define DDMD_ANGLE_POTENTIAL_H

#include <util/param/ParamComposite.h>        // base class
#include <util/boundary/Boundary.h>           // typedef

#include <iostream>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   using namespace Util;

   class Simulation;
   template <int N> class GroupStorage;

   /**
   * An AnglePotential calculates angle forces and energies for a parent Simulation.
   *
   * All operations in this class are local (no MPI).
   */
   class AnglePotential  : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      AnglePotential(Simulation& simulation);

      /**
      * Default constructor (for unit testing).
      */
      AnglePotential();

      /**
      * Destructor.
      */
      ~AnglePotential();

      /**
      * Associate with related objects.
      *
      * Call iff object instantiated with default constructor.
      *
      * \param boundary associated Boundary object.
      * \param storage  associated angle storage object.
      */
      void associate(Boundary& boundary, GroupStorage<3>& storage);

      /// \name Interaction interface
      //@{

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNAngleType(int nAngleType) = 0;
  
      /**
      * Returns potential energy for one angle.
      *
      * \param cosTheta  cosine of the bend angle.
      * \param type      type of bend angle.
      */
      double energy(double cosTheta, int type) const;
 
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
      void force(const Vector& R1, const Vector& R2,
                       Vector& F1, Vector& F2, int type) const;

      #if 0
      /**
      * Return pair interaction class name (e.g., "CosineAngle").
      */
      virtual std::string interactionClassName() const = 0;
      #endif

      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add the angle forces for all atoms.
      */
      virtual void addForces() = 0;

      /**
      * Add pair forces to atom forces, and compute energy.
      */
      virtual void addForces(double& energy) = 0;

      /**
      * Calculate total angle potential.
      *
      * Must be call on all processors (MPI reduce operation).
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator) = 0;
      #else
      virtual void computeEnergy() = 0;
      #endif

      /**
      * Get total angle potential, computed previously by computeEnergy().
      *
      * Call only on master processor, result otherwise invalid.
      */
      virtual double energy() = 0;

      #if 0
      /**
      * Compute total angle pressure.
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const;

      /**
      * Compute x, y, z angle pressure components.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const;

      /**
      * Compute angle stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const;

      //@}
      #endif

   protected:

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated GroupStorage<3> object.
      GroupStorage<3>* storagePtr_;

   };

}
#endif
