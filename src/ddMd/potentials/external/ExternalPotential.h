#ifndef DDMD_EXTERNAL_POTENTIAL_H
#define DDMD_EXTERNAL_POTENTIAL_H

#include <util/param/ParamComposite.h>  // base class
#include <util/boundary/Boundary.h>     // member (typedef)
#include <util/global.h>

#include <iostream>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   class System;
   class AtomStorage;
   class Domain;

   using namespace Util;

   /**
   * An ExternalPotential calculates exteranl forces for a parent System.
   */
   class ExternalPotential : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      ExternalPotential(System& system);

      /**
      * Constructor (for unit testing).
      *
      * \param boundary    associated Boundary object.
      * \param domain      associated Domain object.
      * \param storage     associated AtomStorage object.
      */
      ExternalPotential(Boundary& boundary, Domain& domain, 
                        AtomStorage& storage);

      /**
      * Destructor.
      */
      virtual ~ExternalPotential();

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNAtomType(int nAtomType) = 0;
  
      //@{ External Interaction Interface

      /**
      * Returns external potential energy of a single particle. 
      *
      * \param position atomic position Vector
      * \param i        atom type.
      * \return external potential energy
      */
      virtual double energy(const Vector& position, int i) const = 0;

      /**
      * Returns force caused by the external potential.
      *
      * \param position  atom position
      * \param type      atom type id
      * \param force     force on the atom (on output)
      */
      virtual void getForce(const Vector& position, int type, 
                            Vector& force) const = 0;

      #if 0
      /**
      * Return name of interaction class (e.g., "TanhCosineExternal").
      */
      virtual std::string interactionClassName() const = 0;
      #endif
   
      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add external forces all local atoms.
      */
      virtual void addForces() = 0;

      /**
      * Add pair forces to atom forces, and compute energy.
      */
      virtual void addForces(double& energy) = 0;

      /**
      * Calculate total pair potential on this processor
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator) = 0;
      #else
      virtual void computeEnergy();
      #endif

      /**
      * Return total pair potential on all processors.
      */
      virtual double energy() = 0;

   protected:

      /**
      * Get the parent System by reference.
      */
      System& system();

      /**
      * Get the Boundary by reference.
      */
      Boundary& boundary();

      /**
      * Get the Domain by reference.
      */
      Domain& domain();

      /**
      * Get the AtomStorage by reference.
      */
      AtomStorage& storage();

   private:

      // Pointer to parent System object.
      System* systemPtr_;

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated Domain object.
      Domain* domainPtr_;

      // Pointer to associated AtomStorage object.
      AtomStorage* storagePtr_;

   };

   inline System& ExternalPotential::system() 
   {  return *systemPtr_; }

   inline Boundary& ExternalPotential::boundary() 
   {  return *boundaryPtr_; }

   inline Domain& ExternalPotential::domain()
   {  return *domainPtr_; }

   inline AtomStorage& ExternalPotential::storage()
   {  return *storagePtr_; }

}
#endif
