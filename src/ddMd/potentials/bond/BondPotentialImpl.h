#ifndef DDMD_BOND_POTENTIAL_IMPL_H
#define DDMD_BOND_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BondPotential.h" // base class
#include <util/global.h>

namespace Util
{
   class Vector;
   class Tensor;
}

namespace DdMd
{

   using namespace Util;

   class System;
   template <int N> class GroupStorage;

   /**
   * Implementation template for a BondPotential.
   *
   * \ingroup Bond_Module
   */
   template <class Interaction>
   class BondPotentialImpl : public BondPotential
   {

   public:

      /** 
      * Constructor.
      */
      BondPotentialImpl(System& system);

      /** 
      * Constructor.
      */
      BondPotentialImpl(Boundary& boundary, GroupStorage<2>& storage);

      /** 
      * Destructor.
      */
      virtual ~BondPotentialImpl();

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNBondType(int nAtomType);
  
      /**
      * Read potential energy parameters.
      * 
      * This method reads the bond potential Interaction parameter
      * block. Before calling Interaction::readParam(), it passes
      * simulation().nBondType() to Interaction::setNAtomType().
      */
      virtual void readParam(std::istream& in);

      /// \name Interaction interface
      //@{

      /**
      * Return pair energy for a single pair.
      */
      virtual double energy(double rsq, int bondTypeId) const;

      /**
      * Return force / separation for a single pair.
      */
      virtual double forceOverR(double rsq, int bondTypeId) const;

      /**
      * Return force / separation for a single pair.
      */
      virtual 
      double randomBondLength(Random* random, double beta, int bondTypeId) 
             const;

      #if 0
      /**
      * Return pair interaction class name (e.g., "HarmonicBond").
      */
      virtual std::string interactionClassName() const;
      #endif

      /**
      * Return bond interaction by const reference.
      */
      const Interaction& interaction() const;

      /**
      * Return bond interaction by reference.
      */
      Interaction& interaction();

      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add the bond forces for all atoms.
      */
      virtual void addForces();

      /**
      * Add the bond forces for all atoms and compute energy.
      */
      virtual void addForces(double& energy);

      /**
      * Compute the total bond energy for all processors
      * 
      * Call on all processors (MPI reduce operation).
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator);
      #else
      virtual void computeEnergy();
      #endif

      /**
      * Get the total bond energy, computed previously by computeEnergy().
      *
      * Call only on master. 
      */
      virtual double energy();

      #if 0
      /**
      * Compute total bond pressure.
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const;

      /**
      * Compute x, y, z bond pressure components.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const;

      /**
      * Compute bond stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const;

      //@}
      #endif

   private:

      /**
      * Total bond energy on all processors.
      */
      double energy_;
 
      /**
      * Pointer to Interaction (evaluates the bond potential function).
      */ 
      Interaction* interactionPtr_;

      /*
      * Compute forces and/or energy.
      *
      * Return energy if energy is computed.
      */
      double addForces(bool needForce, bool needEnergy);

      #if 0 
      template <typename T>
      void computeStressImpl(T& stress) const;
      #endif

   };

}

#include "BondPotential.h"
#include <ddMd/system/System.h>
//#include <mcMd/simulation/stress.h>
#include <ddMd/storage/GroupStorage.h>
#include <ddMd/storage/GroupIterator.h>
#include <util/boundary/Boundary.h>
#include <util/space/Vector.h>
#include <util/global.h>

#include <util/space/Dimension.h>
#include <util/space/Vector.h>
//#include <util/space/Tensor.h>
//#include <util/accumulators/setToZero.h>

#include <fstream>

namespace DdMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   template <class Interaction>
   BondPotentialImpl<Interaction>::BondPotentialImpl(System& system)
    : BondPotential(system),
      interactionPtr_(0)
   {  interactionPtr_ = new Interaction(); }
 
   /* 
   * Default constructor.
   */
   template <class Interaction>
   BondPotentialImpl<Interaction>::BondPotentialImpl(Boundary& boundary,
                                   GroupStorage<2>& storage)
    : BondPotential(boundary, storage),
      interactionPtr_(0)
   {  interactionPtr_ = new Interaction(); }
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   BondPotentialImpl<Interaction>::~BondPotentialImpl() 
   {
      if (interactionPtr_) {
         delete interactionPtr_;
         interactionPtr_ = 0;
      }
   }

   /*
   * Set the maximum number of atom types.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::setNBondType(int nBondType)
   {  interaction().setNBondType(nBondType); }

   /*
   * Read bond interaction parameters.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::readParam(std::istream& in)
   {
      readBegin(in,"BondPotential");
      bool nextIndent = false;
      readParamComposite(in, interaction(), nextIndent);
      readEnd(in);
   }

   /*
   * Return bond energy for a single pair.
   */
   template <class Interaction>
   double BondPotentialImpl<Interaction>::energy(double rsq, int bondTypeId) 
      const
   {  return interaction().energy(rsq, bondTypeId); }

   /*
   * Return force / separation for a single bonded pair.
   */
   template <class Interaction>
   double BondPotentialImpl<Interaction>::forceOverR(double rsq, int bondTypeId)
      const
   {  return interaction().forceOverR(rsq, bondTypeId); }

   /*
   * Return random bond length.
   */
   template <class Interaction> double 
   BondPotentialImpl<Interaction>::
      randomBondLength(Random* random, double beta, int bondTypeId) const
   {  return interaction().randomBondLength(random, beta, bondTypeId); }

   #if 0
   /*
   * Return bond potential interaction class name.
   */
   template <class Interaction>
   std::string BondPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }
   #endif

   /**
   * Get Interaction by reference.
   */
   template <class Interaction>
   inline Interaction& BondPotentialImpl<Interaction>::interaction()
   {  return *interactionPtr_; }

   /**
   * Get Interaction by const reference.
   */
   template <class Interaction>
   inline const Interaction& BondPotentialImpl<Interaction>::interaction() const
   {  return *interactionPtr_; }

   /*
   * Increment atomic forces, without calculating energy.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::addForces()
   {  addForces(true, false);  }

   /*
   * Increment atomic forces and compute pair energy for this processor.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::addForces(double& energy)
   {  energy = addForces(true, true);  }

   /*
   * Compute total bond energy on all processors.
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void 
   BondPotentialImpl<Interaction>::computeEnergy(MPI::Intracomm& communicator)
   #else
   void BondPotentialImpl<Interaction>::computeEnergy()
   #endif
   { 
      double localEnergy = 0; 
      localEnergy = addForces(false, true); 
      #ifdef UTIL_MPI
      communicator.Reduce(&localEnergy, &energy_, 1, 
                          MPI::DOUBLE, MPI::SUM, 0);
      #else
      energy_ = localEnergy;
      #endif
   }

   /*
   * Return total pair energy from all processors.
   */
   template <class Interaction>
   double BondPotentialImpl<Interaction>::energy()
   {  return energy_; } 

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction> double 
   BondPotentialImpl<Interaction>::addForces(bool needForce, bool needEnergy)
   {
      // Preconditions
      //if (!storagePtr_->isInitialized()) {
      //   UTIL_THROW("AtomStorage must be initialized");
      //}

      Vector f;
      double rsq;
      double bondEnergy;
      double energy = 0.0;
      GroupIterator<2> iter;
      Atom* atom0Ptr;
      Atom* atom1Ptr;
      int type;
      int isLocal0, isLocal1;

      storagePtr_->begin(iter);
      for ( ; !iter.atEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         // Set f = r0 - r1, minimum image separation between atoms
         rsq = boundaryPtr_->distanceSq(atom0Ptr->position(), 
                                        atom1Ptr->position(), f);
         if (needEnergy) {
            bondEnergy = interactionPtr_->energy(rsq, type);
            if (isLocal0 && isLocal1) {
               energy += bondEnergy;
            } else {
               energy += 0.5*bondEnergy;
            }
         }
         if (needForce) {
            // Set force = (r0-r1)*(forceOverR)
            f *= interactionPtr_->forceOverR(rsq, type);
            if (isLocal0) {
               atom0Ptr->force() += f;
            }
            if (isLocal1) {
               atom1Ptr->force() -= f;
            }
         }
      }
      return energy;
   }

   #if 0
   /* 
   * Compute total stress (generic).
   */
   template <class Interaction>
   template <typename T>
   void BondPotentialImpl<Interaction>::computeStressImpl(T& stress) 
      const
   {
      Vector dr, force;
      double rsq;
      const Atom* atom0Ptr;
      const Atom* atom1Ptr;

      setToZero(stress);

      // Iterate over all bonds in System
      System::ConstMoleculeIterator molIter;
      Molecule::ConstBondIterator bondIter;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (simulation().species(iSpec).nBond() > 0) {
            for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               molIter->begin(bondIter); 
               for ( ; bondIter.notEnd(); ++bondIter) {
                  atom0Ptr = &(bondIter->atom(0));
                  atom1Ptr = &(bondIter->atom(1));
                  rsq = boundary().distanceSq(atom0Ptr->position(), 
                                              atom1Ptr->position(), dr);
                  force = dr;
                  force *= interaction().forceOverR(rsq, 
                                                  bondIter->typeId());
                  incrementPairStress(force, dr, stress);
               }
            }
         }
      }

      stress /= boundary().volume();
      normalizeStress(stress);
   }

   template <class Interaction>
   void BondPotentialImpl<Interaction>::computeStress(double& stress) const
   {  computeStressImpl(stress); }

   template <class Interaction>
   void BondPotentialImpl<Interaction>::computeStress(Util::Vector& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Interaction>
   void BondPotentialImpl<Interaction>::computeStress(Util::Tensor& stress) 
        const
   {  computeStressImpl(stress); }
   #endif

}
#endif
