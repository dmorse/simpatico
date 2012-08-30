#ifndef DDMD_EXTERNAL_POTENTIAL_H
#define DDMD_EXTERNAL_POTENTIAL_H

#include <util/param/ParamComposite.h>  // base class
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
   class ExternalPotential : public ParamComposite
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
      virtual void computeForces() = 0;

      /**
      * Add pair forces to atom forces, and compute energy.
      */
      virtual void computeForces(double& energy) = 0;

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
