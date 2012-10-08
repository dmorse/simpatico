#ifndef DDMD_EXTERNAL_POTENTIAL_H
#define DDMD_EXTERNAL_POTENTIAL_H

#include <ddMd/potentials/Potential.h>  // base class
#include <util/boundary/Boundary.h>     // member (typedef)
#include <util/global.h>

#include <iostream>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   class Simulation;
   class AtomStorage;

   using namespace Util;

   /**
   * An ExternalPotential calculates exteranl forces for a parent Simulation.
   *
   * \ingroup DdMd_External_Module
   */
   class ExternalPotential : public Potential
   {

   public:

      /**
      * Constructor.
      */
      ExternalPotential(Simulation& simulation);

      /**
      * Default constructor (for unit testing).
      */
      ExternalPotential();

      /**
      * Associate with related objects.
      *
      * Call iff object instantiated with default constructor.
      *
      * \param boundary associated Boundary object.
      * \param storage  associated AtomStorage object.
      */
      virtual void associate(Boundary& boundary, AtomStorage& storage);

      /**
      * Destructor.
      */
      virtual ~ExternalPotential();

      /// \name Interaction interface
      //@{

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNAtomType(int nAtomType) = 0;
  
      /**
      * Returns external potential energy of a single particle. 
      *
      * \param position atomic position Vector
      * \param i        atom type.
      * \return external potential energy
      */
      virtual double externalEnergy(const Vector& position, int i) const = 0;

      /**
      * Returns force caused by the external potential.
      *
      * \param position  atom position
      * \param type      atom type id
      * \param force     force on the atom (on output)
      */
      virtual void getExternalForce(const Vector& position, int type, 
                                    Vector& force) const = 0;

      /**
      * Return name of interaction class (e.g., "TanhCosineExternal").
      */
      virtual std::string interactionClassName() const = 0;
   
      //@}

   protected:

      /**
      * Get the parent Simulation by reference.
      */
      Simulation& simulation();

      /**
      * Get the Boundary by reference.
      */
      Boundary& boundary();

      /**
      * Get the AtomStorage by reference.
      */
      AtomStorage& storage();

   private:

      // Pointer to parent Simulation object.
      Simulation* simulationPtr_;

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated AtomStorage object.
      AtomStorage* storagePtr_;

   };

   inline Simulation& ExternalPotential::simulation() 
   {  return *simulationPtr_; }

   inline Boundary& ExternalPotential::boundary() 
   {  return *boundaryPtr_; }

   inline AtomStorage& ExternalPotential::storage()
   {  return *storagePtr_; }

}
#endif
