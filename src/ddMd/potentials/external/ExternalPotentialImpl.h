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
      * Returns external potential energy of a single particle. 
      *
      * \param position atomic position Vector
      * \param typeId   atom type index
      * \return external potential energy
      */
      virtual double energy(const Vector& position, int typeId) const;

      /**
      * Returns force caused by the external potential.
      *
      * \param position  atom position
      * \param typeId    atom type index
      * \param force     force on the atom (on output)
      */
      virtual void getForce(const Vector& position, int typeId, 
                            Vector& force) const; 

      #if 0
      /**
      * Return external interaction class name (e.g., "TanhCosineExternal").
      */
      virtual std::string interactionClassName() const;
      #endif

      /**
      * Return external interaction by reference.
      */
      Interaction& interaction();

      /**
      * Return external interaction by const reference.
      */
      const Interaction& interaction() const;

      //@}
      /// \name Total Energy, Forces and Stress 
      //@{

      /**
      * Calculate external forces for all atoms in this Simulation.
      */
      virtual void computeForces();

      /**
      * Add external forces to atom forces, and compute energy.
      *
      * \param energy on output, contains energy for this processor.
      */
      virtual void computeForces(double& energy);

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

      /**
      * Get the total energy calculated previously by computeEnergy().
      *
      * Call only on master. 
      */
      virtual double energy();

      #if 0 
      /**
      * Compute total nonbonded pressure
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const;

      /**
      * Compute x, y, z nonbonded pressures.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const;

      /**
      * Compute nonbonded stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const;
      #endif

      //@}

   private:

      /**
      * Total external energy on all processors (value valid only on master).
      */
      double energy_;
 
      /**
      * Pointer to external interaction object.
      */ 
      Interaction* interactionPtr_;

      /**
      * Calculate external forces and/or energy.
      */
      double computeForces(bool needForce, bool needEnergy);

      #if 0 
      template <typename T>
      void computeStressImpl(T& stress) const;
      #endif

   };

}

#include <ddMd/simulation/Simulation.h>
//#include <ddMd/simulation/stress.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>

#include <util/space/Dimension.h>
#include <util/space/Vector.h>
//#include <util/space/Tensor.h>
//#include <util/accumulators/setToZero.h>
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
   ExternalPotentialImpl<Interaction>::energy(const Vector& position, int typeId) 
      const
   {  return interaction().energy(position, typeId); }

   /*
   * Return external force on an atom.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::getForce(const Vector& position, 
                                               int typeId, Vector& force) const
   { interaction().getForce(position, typeId, force); }

   #if 0
   /*
   * Return external interaction class name.
   */
   template <class Interaction>
   std::string ExternalPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }
   #endif

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
   * Increment atomic forces, and compute total external energy.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::computeForces(double& energy)
   {  
      energy = computeForces(true, true); 
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
      communicator.Reduce(&localEnergy, &energy_, 1, 
                          MPI::DOUBLE, MPI::SUM, 0);
      #else
      energy_ = localEnergy;
      #endif
   }

   /*
   * Return total external energy from all processors.
   */
   template <class Interaction>
   double ExternalPotentialImpl<Interaction>::energy()
   {  return energy_; } 

   /*
   * Increment atomic forces and/or external energy (private).
   */
   template <class Interaction>
   double 
   ExternalPotentialImpl<Interaction>::computeForces(bool needForce, 
                                                 bool needEnergy)
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
