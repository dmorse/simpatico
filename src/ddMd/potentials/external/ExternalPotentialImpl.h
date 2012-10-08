#ifndef DDMD_EXTERNAL_POTENTIAL_IMPL_H
#define DDMD_EXTERNAL_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/potentials/external/ExternalPotential.h>
#include <util/global.h>

namespace Util {
   class Vector;
   class Tensor;
}

namespace DdMd
{

   using namespace Util;

   class Simulation;
   class AtomStorage;

   /**
   * Implementation template for a ExternalPotential.
   *
   * \ingroup DdMd_External_Module
   */
   template <class Interaction>
   class ExternalPotentialImpl : public ExternalPotential
   {

   public:

      /** 
      * Constructor.
      */
      ExternalPotentialImpl(Simulation& simulation);

      /** 
      * Default constructor (for unit testing).
      */
      ExternalPotentialImpl();

      /** 
      * Destructor.
      */
      virtual ~ExternalPotentialImpl();

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
      * Read external potential interaction.
      * 
      * \param in input parameter stream.
      */
      virtual void readParameters(std::istream& in);

      /// \name Interaction interface
      //@{

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNAtomType(int nAtomType);

      /**
      * Returns external potential energy of a single atom.
      *
      * \param position atomic position Vector
      * \param typeId   atom type index
      * \return external potential energy
      */
      virtual double externalEnergy(const Vector& position, int typeId) const;

      /**
      * Computes external force on one atom.
      *
      * \param position  atom position (input)
      * \param typeId    atom type index (input)
      * \param force     force on the atom (output)
      */
      virtual void getExternalForce(const Vector& position, int typeId, 
                                    Vector& force) const; 

      /**
      * Return external interaction class name (e.g., "TanhCosineExternal").
      */
      virtual std::string interactionClassName() const;

      /**
      * Return external interaction by reference.
      */
      Interaction& interaction();

      /**
      * Return external interaction by const reference.
      */
      const Interaction& interaction() const;

      //@}
      /// \name Total Energy and Forces 
      //@{

      /**
      * Calculate external forces for all atoms in this Simulation.
      */
      virtual void computeForces();

      /**
      * Compute the total external energy for all processors
      * 
      * Call on all processors.
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator);
      #else
      virtual void computeEnergy();
      #endif

      //@}

   private:

      /**
      * Pointer to external interaction object.
      */ 
      Interaction* interactionPtr_;

      /**
      * Calculate external forces and/or energy.
      */
      double computeForces(bool needForce, bool needEnergy);

   };

}

#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>

#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/global.h>

#include <fstream>

namespace DdMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   template <class Interaction>
   ExternalPotentialImpl<Interaction>::ExternalPotentialImpl(Simulation& simulation)
    : ExternalPotential(simulation),
      interactionPtr_(0)
   {  
      interactionPtr_ = new Interaction;
      //interaction().setNAtomType(simulation.nAtomType());
      interactionPtr_->setBoundary(simulation.boundary());
   }
 
   /* 
   * Default constructor.
   */
   template <class Interaction>
   ExternalPotentialImpl<Interaction>::ExternalPotentialImpl()
    : ExternalPotential(),
      interactionPtr_(0)
   {  interactionPtr_ = new Interaction; }
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   ExternalPotentialImpl<Interaction>::~ExternalPotentialImpl() 
   {
      if (interactionPtr_) {
         delete interactionPtr_;
         interactionPtr_ = 0;
      }
   }

   /*
   * Associate with related objects. (for unit testing).
   *
   * Required iff instantiated with default constructor.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>
        ::associate(Boundary& boundary, AtomStorage& storage)
   {
      ExternalPotential::associate(boundary, storage);
      interactionPtr_->setBoundary(boundary);
   } 

   /* 
   * Read parameters from file.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::readParameters(std::istream &in) 
   {
      bool nextIndent = false; // Do not indent interaction block. 
      addParamComposite(interaction(), nextIndent);
      interaction().readParameters(in);
   }
  
   /**
   * Set the maximum number of atom types.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::setNAtomType(int nAtomType)
   {  interaction().setNAtomType(nAtomType); }

   /*
   * Return external energy of an atom.
   */
   template <class Interaction>
   double 
   ExternalPotentialImpl<Interaction>::externalEnergy(const Vector& position, 
                                                      int typeId) const
   {  return interaction().energy(position, typeId); }

   /*
   * Return external force on an atom.
   */
   template <class Interaction>
   void 
    ExternalPotentialImpl<Interaction>::getExternalForce(const Vector& position, 
                                                int typeId, Vector& force) const
   { interaction().getForce(position, typeId, force); }

   /*
   * Return external interaction class name.
   */
   template <class Interaction>
   std::string ExternalPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }

   /**
   * Return underlying interaction object by const reference.
   */
   template <class Interaction>
   const Interaction& ExternalPotentialImpl<Interaction>::interaction() const
   {  return *interactionPtr_; }

   /**
   * Return underlying interaction object by reference.
   */
   template <class Interaction>
   Interaction& ExternalPotentialImpl<Interaction>::interaction()
   {  return *interactionPtr_; }

   /*
   * Increment atomic forces, without calculating energy.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::computeForces()
   {  
      computeForces(true, false); 
   }

   /*
   * Compute total external energy on all processors.
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void 
   ExternalPotentialImpl<Interaction>::computeEnergy(MPI::Intracomm& communicator)
   #else
   void ExternalPotentialImpl<Interaction>::computeEnergy()
   #endif
   { 
      double localEnergy = 0; 
      localEnergy = computeForces(false, true); 

      #ifdef UTIL_MPI
      reduceEnergy(localEnergy, communicator);
      #else
      setEnergy(localEnergy);
      #endif
   }

   /*
   * Increment atomic forces and/or external energy (private).
   */
   template <class Interaction>
   double 
   ExternalPotentialImpl<Interaction>::computeForces(bool needForce, bool needEnergy)
   {
      // Preconditions
      //if (!storage().isInitialized()) {
      //   UTIL_THROW("AtomStorage must be initialized");
      //}

      Vector f;
      double energy = 0.0;
      AtomIterator iter;
      int type;

      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         if (needEnergy) {
            energy += 
               interactionPtr_->energy(iter->position(), type);
         }
         if (needForce) {
            interactionPtr_->getForce(iter->position(), type, f);
            iter->force() += f;
         }
      }
      return energy;
   }

}
#endif
