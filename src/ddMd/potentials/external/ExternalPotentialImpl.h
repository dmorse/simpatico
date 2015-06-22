#ifndef DDMD_EXTERNAL_POTENTIAL_IMPL_H
#define DDMD_EXTERNAL_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
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
      * Set the maximum number of atom types.
      *
      * Call iff object instantiated with default constructor.
      *
      * \param nAtomType maximum number of atomTypes (max index+1)
      */
      virtual void setNAtomType(int nAtomType);

      /**
      * Read external potential interaction.
      * 
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from archive and initialize, for restart.
      * 
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive, for restart.
      *
      * Call only on master processor.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /// \name Interaction interface
      //@{

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
      * Modify a pair interaction parameter, identified by a string.
      *
      * \param name   parameter name
      * \param value  new value of parameter
      */
      void set(std::string name, double value)
      {  interactionPtr_->set(name, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      */
      double get(std::string name) const
      {  return interactionPtr_->get(name); }

      /**
      * Return external interaction class name (e.g., "LamellarOrderingExternal").
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
      * This is initialized to false, and set true in (read|load)Parameters.
      */
      bool isInitialized_;

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
      interactionPtr_(0),
      isInitialized_(false)
   {  
      interactionPtr_ = new Interaction;
      interaction().setBoundary(simulation.boundary());
      setNAtomType(simulation.nAtomType());
   }
 
   /* 
   * Default constructor.
   */
   template <class Interaction>
   ExternalPotentialImpl<Interaction>::ExternalPotentialImpl()
    : ExternalPotential(),
      interactionPtr_(0),
      isInitialized_(false)
   {  interactionPtr_ = new Interaction; }
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   ExternalPotentialImpl<Interaction>::~ExternalPotentialImpl() 
   {
      if (interactionPtr_) {
         delete interactionPtr_;
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
      interaction().setBoundary(boundary);
   } 

   /**
   * Set the maximum number of atom types (for unit testing).
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::setNAtomType(int nAtomType)
   {  interaction().setNAtomType(nAtomType); }

   /* 
   * Read parameters from file.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::readParameters(std::istream &in) 
   {
      UTIL_CHECK(!isInitialized_);
      bool nextIndent = false; // Do not indent interaction block
      addParamComposite(interaction(), nextIndent);
      interaction().readParameters(in);
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::loadParameters(Serializable::IArchive &ar)
   {
      UTIL_CHECK(!isInitialized_);
      bool nextIndent = false; // Do not indent interaction block
      addParamComposite(interaction(), nextIndent);
      interaction().loadParameters(ar);
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive (call only on master processor).
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::save(Serializable::OArchive &ar)
   {  interaction().save(ar); }
 
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
   {
      assert(interactionPtr_);  
      return *interactionPtr_; 
   }

   /**
   * Return underlying interaction object by reference.
   */
   template <class Interaction>
   Interaction& ExternalPotentialImpl<Interaction>::interaction()
   {  
      assert(interactionPtr_);  
      return *interactionPtr_; 
   }

   /*
   * Increment atomic forces, without calculating energy.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::computeForces()
   {  computeForces(true, false); }

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
