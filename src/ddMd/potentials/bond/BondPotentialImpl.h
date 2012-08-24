#ifndef DDMD_BOND_POTENTIAL_IMPL_H
#define DDMD_BOND_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BondPotential.h" // base class
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   class Simulation;
   template <int N> class GroupStorage;

   /**
   * Implementation template for a BondPotential.
   *
   * \ingroup DdMd_Bond_Module
   */
   template <class Interaction>
   class BondPotentialImpl : public BondPotential
   {

   public:

      /** 
      * Constructor.
      */
      BondPotentialImpl(Simulation& simulation);

      /** 
      * Default constructor.
      */
      BondPotentialImpl();

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
      virtual double bondEnergy(double rsq, int bondTypeId) const;

      /**
      * Return force / separation for a single pair.
      */
      virtual double bondForceOverR(double rsq, int bondTypeId) const;

      /**
      * Return force / separation for a single pair.
      */
      virtual 
      double randomBondLength(Random* random, double beta, int bondTypeId) 
             const;

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
      * Compute the covalent bond stress.
      * 
      * Call on all processors.
      */
      #ifdef UTIL_MPI
      virtual void computeStress(MPI::Intracomm& communicator);
      #else
      virtual void computeStress();
      #endif

      //@}

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

   };

}

#include "BondPotential.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/GroupStorage.h>
#include <ddMd/storage/GroupIterator.h>
#include <util/boundary/Boundary.h>
#include <util/space/Vector.h>
#include <util/global.h>

#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
//#include <util/accumulators/setToZero.h>

#include <fstream>

namespace DdMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   template <class Interaction>
   BondPotentialImpl<Interaction>::BondPotentialImpl(Simulation& simulation)
    : BondPotential(simulation),
      interactionPtr_(0)
   {  interactionPtr_ = new Interaction(); }
 
   /* 
   * Default constructor.
   */
   template <class Interaction>
   BondPotentialImpl<Interaction>::BondPotentialImpl()
    : BondPotential(),
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
   double 
   BondPotentialImpl<Interaction>::bondEnergy(double rsq, int bondTypeId) const
   {  return interaction().energy(rsq, bondTypeId); }

   /*
   * Return force / separation for a single bonded pair.
   */
   template <class Interaction>
   double 
   BondPotentialImpl<Interaction>::bondForceOverR(double rsq, int bondTypeId)
   const
   {  return interaction().forceOverR(rsq, bondTypeId); }

   /*
   * Return random bond length.
   */
   template <class Interaction> double 
   BondPotentialImpl<Interaction>::
      randomBondLength(Random* random, double beta, int bondTypeId) const
   {  return interaction().randomBondLength(random, beta, bondTypeId); }

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

      // Do nothing and return if energy is already set.
      if (isEnergySet()) return;
 
      double localEnergy = 0.0; 
      localEnergy = addForces(false, true); 
      #ifdef UTIL_MPI
      double totalEnergy = 0.0; 
      communicator.Reduce(&localEnergy, &totalEnergy, 1, 
                          MPI::DOUBLE, MPI::SUM, 0);
      if (communicator.Get_rank() != 0) {
         totalEnergy = 0.0;
      }
      setEnergy(totalEnergy);
      #else
      setEnergy(localEnergy);
      #endif
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction> double 
   BondPotentialImpl<Interaction>::addForces(bool needForce, bool needEnergy)
   {
      // Preconditions
      //if (!storage().isInitialized()) {
      //   UTIL_THROW("BondStorage must be initialized");
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

      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         // Set f = r0 - r1, minimum image separation between atoms
         rsq = boundary().distanceSq(atom0Ptr->position(), 
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

   /*
   * Compute total pair stress (Call on all processors).
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void BondPotentialImpl<Interaction>::computeStress(MPI::Intracomm& communicator)
   #else
   void BondPotentialImpl<Interaction>::computeStress()
   #endif
   {
      // Do nothing and return if stress is already set.
      if (isStressSet()) return;
 
      Tensor localStress;
      Vector dr;
      Vector f;
      double rsq, forceOverR;
      GroupIterator<2> iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int    type;
      int    isLocal0, isLocal1;

      localStress.zero();

      // Iterate over bonds
      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         rsq = boundary().distanceSq(atom0Ptr->position(), 
                                        atom1Ptr->position(), dr);
         f = dr;
         assert(isLocal0 || isLocal1);
         if (isLocal0 && isLocal1) {
            forceOverR = interactionPtr_->forceOverR(rsq, type);
         } else {
            forceOverR = 0.5*interactionPtr_->forceOverR(rsq, type);
         }
         f *= forceOverR;
         incrementPairStress(f, dr, localStress);
      }

      // if (reverseUpdateFlag()) { } else { } 

      // Normalize by volume 
      localStress /= boundary().volume();

      #ifdef UTIL_MPI
      // Reduce results from all processors
      Tensor totalStress;
      communicator.Reduce(&localStress(0, 0), &totalStress(0, 0), 
                          Dimension*Dimension, MPI::DOUBLE, MPI::SUM, 0);
      if (communicator.Get_rank() != 0) {
         totalStress.zero();
      }
      setStress(totalStress);
      #else
      setStress(localStress);
      #endif

   }

}
#endif
