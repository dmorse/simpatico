#ifndef DDMD_DIHEDRAL_POTENTIAL_H
#define DDMD_DIHEDRAL_POTENTIAL_H

#include <ddMd/potentials/Potential.h>       // base class
#include <util/boundary/Boundary.h>          // typedef

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
   * \ingroup DdMd_Dihedral_Module
   */
   class DihedralPotential  : public Potential
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
      virtual double 
      dihedralEnergy(const Vector& R1, const Vector& R2, const Vector& R3,
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
      virtual void 
      dihedralForce(const Vector& R1, const Vector& R2, const Vector& R3,
                    Vector& F1, Vector& F2, Vector& F3, int type) const = 0;

      /**
      * Modify a dihedral interaction parameter, identified by a string.
      *
      * \param name  parameter variable name
      * \param type  dihedral type index
      * \param value  new value of parameter
      */
      virtual void set(std::string name, int type, double value) = 0;

      /**
      * Get a dihedral parameter value, identified by a string.
      *
      * \param name  parameter variable name
      * \param type  dihedral type index
      */
      virtual double get(std::string name, int type) const = 0;

      /**
      * Return pair interaction class name (e.g., "HarmonicDihedral").
      */
      virtual std::string interactionClassName() const = 0;

      //@}

   protected:

      /**
      *  Return boundary by reference.   
      */
      Boundary& boundary() const;

      /**
      *  Return bond storage by reference.   
      */
      GroupStorage<4>& storage() const;

   private:

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated GroupStorage<4> object.
      GroupStorage<4>* storagePtr_;

   };

   // Get boundary by reference.
   inline Boundary& DihedralPotential::boundary() const
   { return *boundaryPtr_; }

   // Get bond storage by reference.
   inline GroupStorage<4>& DihedralPotential::storage() const
   { return *storagePtr_; }


}
#endif
