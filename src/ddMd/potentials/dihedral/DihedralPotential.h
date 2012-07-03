#ifndef DDMD_DIHEDRAL_POTENTIAL_H
#define DDMD_DIHEDRAL_POTENTIAL_H

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
   * An DihedralPotential calculates dihedral forces and energies for a parent Simulation.
   *
   * All operations in this class are local (no MPI).
   */
   class DihedralPotential  : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      DihedralPotential(Simulation& simulation);

      /**
      * Default constructor (for unit testing).
      */
      DihedralPotential();

      /**
      * Destructor.
      */
      ~DihedralPotential();

      /**
      * Associate with related objects.
      *
      * Call iff object instantiated with default constructor.
      *
      * \param boundary associated Boundary object.
      * \param storage  associated dihedral storage object.
      */
      void associate(Boundary& boundary, GroupStorage<4>& storage);

      /// \name Interaction interface
      //@{

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNDihedralType(int nDihedralType) = 0;
  
      /**
      * Returns potential energy for one dihedral.
      *
      *     0   2    3
      *     o   o----o
      *      \ /
      *       o 
      *       1 
      *
      * \param R1     bond vector from atom 0 to 1
      * \param R2     bond vector from atom 1 to 2
      * \param R3     bond vector from atom 2 to 3
      * \param type   type of dihedral
      */
      virtual
      double energy(const Vector& R1, const Vector& R2, const Vector& R3,
                    int type) const = 0;
 
      /**
      * Returns derivatives of energy with respect to bond vectors forming the
      * dihedral group.
      *
      * \param R1     bond vector from atom 0 to 1
      * \param R2     bond vector from atom 1 to 2
      * \param R3     bond vector from atom 2 to 3
      * \param F1     force along R1 direction (upon return)
      * \param F2     force along R2 direction (upon return)
      * \param F3     force along R2 direction (upon return)
      * \param type   type of dihedral
      */
      virtual
      void force(const Vector& R1, const Vector& R2, const Vector& R3,
                 Vector& F1, Vector& F2, Vector& F3, int type) const = 0;

      #if 0
      /**
      * Return pair interaction class name (e.g., "HarmonicDihedral").
      */
      virtual std::string interactionClassName() const = 0;
      #endif

      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add the dihedral forces for all atoms.
      */
      virtual void addForces() = 0;

      /**
      * Add pair forces to atom forces, and compute energy.
      */
      virtual void addForces(double& energy) = 0;

      /**
      * Calculate total dihedral potential.
      *
      * Must be call on all processors (MPI reduce operation).
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator) = 0;
      #else
      virtual void computeEnergy() = 0;
      #endif

      /**
      * Get total dihedral potential, computed previously by computeEnergy().
      *
      * Call only on master processor, result otherwise invalid.
      */
      virtual double energy() = 0;

      #if 0
      /**
      * Compute total dihedral pressure.
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const;

      /**
      * Compute x, y, z dihedral pressure components.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const;

      /**
      * Compute dihedral stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const;

      //@}
      #endif

   protected:

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated GroupStorage<4> object.
      GroupStorage<4>* storagePtr_;

   };

}
#endif
