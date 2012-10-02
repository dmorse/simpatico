#ifndef DDMD_ANGLE_POTENTIAL_H
#define DDMD_ANGLE_POTENTIAL_H

#include <ddMd/potentials/Potential.h>        // base class
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
   * \ingroup DdMd_Angle_Module
   */
   class AnglePotential  : public Potential
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
      * Set the maximum number of angle types.
      */
      virtual void setNAngleType(int nAngleType) = 0;
  
      /**
      * Returns potential energy for one angle.
      *
      * \param cosTheta  cosine of the angle.
      * \param type  type of angle group
      */
      double angleEnergy(double cosTheta, int type) const;
 
      /**
      * Computes forces along two bonds within the angle.
      *
      * \param R1     bond vector from atom 1 to 2.
      * \param R2     bond vector from atom 2 to 3.
      * \param F1     return force along R1 direction.
      * \param F2     return force along R2 direction.
      * \param type   type of angle.
      */
      void angleForce(const Vector& R1, const Vector& R2,
                       Vector& F1, Vector& F2, int type) const;

      /**
      * Modify a angle interaction parameter, identified by a string.
      *
      * \param name  parameter variable name
      * \param type  angle type index
      * \param value  new value of parameter
      */
      virtual void set(std::string name, int type, double value) = 0;

      /**
      * Get a angle parameter value, identified by a string.
      *
      * \param name  parameter variable name
      * \param type  angle type index
      */
      virtual double get(std::string name, int type) const = 0;

      /**
      * Return pair interaction class name (e.g., "CosineAngle").
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
      GroupStorage<3>& storage() const;

   private:

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated GroupStorage<3> object.
      GroupStorage<3>* storagePtr_;

   };

   // Get boundary by reference.
   inline Boundary& AnglePotential::boundary() const
   { return *boundaryPtr_; }

   // Get bond storage by reference.
   inline GroupStorage<3>& AnglePotential::storage() const
   { return *storagePtr_; }

}
#endif
