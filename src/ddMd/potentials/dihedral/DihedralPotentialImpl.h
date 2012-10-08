#ifndef DDMD_DIHEDRAL_POTENTIAL_IMPL_H
#define DDMD_DIHEDRAL_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DihedralPotential.h" // base class
#include <util/global.h>

namespace Util
{
   class Vector;
   class Tensor;
}

namespace DdMd
{

   using namespace Util;

   class Simulation;
   template <int N> class GroupStorage;

   /**
   * Implementation template for a DihedralPotential.
   *
   * \ingroup DdMd_Dihedral_Module
   */
   template <class Interaction>
   class DihedralPotentialImpl : public DihedralPotential
   {

   public:

      /** 
      * Constructor.
      */
      DihedralPotentialImpl(Simulation& simulation);

      /** 
      * Default constructor.
      */
      DihedralPotentialImpl();

      /** 
      * Destructor.
      */
      virtual ~DihedralPotentialImpl();

      /**
      * Read potential energy parameters.
      * 
      * This method reads the dihedral potential Interaction parameter
      * block. Before calling Interaction::readParameters(), it passes
      * simulation().nDihedralType() to Interaction::setNAtomType().
      */
      virtual void readParameters(std::istream& in);

      /// \name Interaction interface
      //@{

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNDihedralType(int nAtomType);
  
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
                     int type) const;
 
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
                    Vector& F1, Vector& F2, Vector& F3, int type) const;

      /**
      * Modify an dihedral parameter, identified by a string.
      *
      * \param name  parameter variable name
      * \param type  type index for dihedral group
      * \param value  new value of parameter
      */
      void set(std::string name, int type, double value)
      {  interactionPtr_->set(name, type, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name  parameter variable name
      * \param type  type index of dihedral group
      */
      double get(std::string name, int type) const
      {  return interactionPtr_->get(name, type); }

      /**
      * Return pair interaction class name (e.g., "CosineDihedral").
      */
      virtual std::string interactionClassName() const;

      /**
      * Return dihedral interaction by const reference.
      */
      const Interaction& interaction() const;

      /**
      * Return dihedral interaction by reference.
      */
      Interaction& interaction();

      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add the dihedral forces for all atoms.
      */
      virtual void computeForces();

      /**
      * Compute the total dihedral energy for all processors
      * 
      * Call on all processors (MPI reduce operation).
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator);
      #else
      virtual void computeEnergy();
      #endif

      /**
      * Compute the covalent dihedral stress.
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
      * Pointer to Interaction (evaluates the dihedral potential function).
      */ 
      Interaction* interactionPtr_;

      /*
      * Compute forces and/or energy.
      *
      * Return energy if energy is computed.
      */
      double computeForces(bool needForce, bool needEnergy);

   };

}

#include "DihedralPotential.h"
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
   DihedralPotentialImpl<Interaction>::DihedralPotentialImpl(Simulation& simulation)
    : DihedralPotential(simulation),
      interactionPtr_(0)
   {  interactionPtr_ = new Interaction(); }
 
   /* 
   * Default constructor.
   */
   template <class Interaction>
   DihedralPotentialImpl<Interaction>::DihedralPotentialImpl()
    : DihedralPotential(),
      interactionPtr_(0)
   {  interactionPtr_ = new Interaction(); }
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   DihedralPotentialImpl<Interaction>::~DihedralPotentialImpl() 
   {
      if (interactionPtr_) {
         delete interactionPtr_;
         interactionPtr_ = 0;
      }
   }

   /*
   * Set the maximum number of dihedral types.
   */
   template <class Interaction>
   void DihedralPotentialImpl<Interaction>::setNDihedralType(int nDihedralType)
   {  interaction().setNDihedralType(nDihedralType); }

   /*
   * Read dihedral interaction parameters.
   */
   template <class Interaction>
   void DihedralPotentialImpl<Interaction>::readParameters(std::istream& in)
   {
      bool nextIndent = false;
      addParamComposite(interaction(), nextIndent);
      interaction().readParameters(in);
   }

   /*
   * Return energy for a single dihedral.
   */
   template <class Interaction>
   inline double DihedralPotentialImpl<Interaction>::
      dihedralEnergy(const Vector& R1, const Vector& R2, const Vector& R3, 
                     int typeId) const
   {  return interaction().energy(R1, R2, R3, typeId); }

   /*
   * Return forces for a single dihedral.
   */
   template <class Interaction>
   inline void DihedralPotentialImpl<Interaction>::
      dihedralForce(const Vector& R1, const Vector& R2, const Vector& R3,
                    Vector& F1, Vector& F2, Vector& F3, int typeId) const
   {  interaction().force(R1, R2, R3, F1, F2, F3, typeId); }

   /*
   * Return dihedral potential interaction class name.
   */
   template <class Interaction>
   std::string DihedralPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }

   /**
   * Get Interaction by reference.
   */
   template <class Interaction>
   inline Interaction& DihedralPotentialImpl<Interaction>::interaction()
   {  return *interactionPtr_; }

   /**
   * Get Interaction by const reference.
   */
   template <class Interaction>
   inline const Interaction& DihedralPotentialImpl<Interaction>::interaction() const
   {  return *interactionPtr_; }

   /*
   * Increment atomic forces, without calculating energy.
   */
   template <class Interaction>
   void DihedralPotentialImpl<Interaction>::computeForces()
   {  computeForces(true, false);  }

   #if 0
   /*
   * Increment atomic forces and compute pair energy for this processor.
   */
   template <class Interaction>
   void DihedralPotentialImpl<Interaction>::computeForces(double& energy)
   {  energy = computeForces(true, true);  }
   #endif

   /*
   * Compute total dihedral energy on all processors.
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void 
   DihedralPotentialImpl<Interaction>::computeEnergy(MPI::Intracomm& communicator)
   #else
   void DihedralPotentialImpl<Interaction>::computeEnergy()
   #endif
   { 

      // Do nothing and return if energy is already set.
      if (isEnergySet()) return;
 
      double localEnergy = 0.0; 
      localEnergy = computeForces(false, true); 

      reduceEnergy(localEnergy, communicator);
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction> double 
   DihedralPotentialImpl<Interaction>::computeForces(bool needForce, bool needEnergy)
   {
      // Preconditions
      //if (!storage().isInitialized()) {
      //   UTIL_THROW("GroupStorage must be initialized");
      //}

      Vector dr1; // R[1] - R[0]
      Vector dr2; // R[2] - R[1]
      Vector dr3; // R[3] - R[2]
      Vector f1, f2, f3;
      double energy = 0.0;
      double dihedralEnergy;
      double fraction;
      GroupIterator<4> iter;
      Atom* atom0Ptr;
      Atom* atom1Ptr;
      Atom* atom2Ptr;
      Atom* atom3Ptr;
      int   type, isLocal0, isLocal1, isLocal2, isLocal3;

      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         atom2Ptr = iter->atomPtr(2);
         atom3Ptr = iter->atomPtr(3);
         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         isLocal2 = !(atom2Ptr->isGhost());
         isLocal3 = !(atom3Ptr->isGhost());
         // Calculate minimimum image separations
         boundary().distanceSq(atom1Ptr->position(),
                                  atom0Ptr->position(), dr1);
         boundary().distanceSq(atom2Ptr->position(),
                                  atom1Ptr->position(), dr2);
         boundary().distanceSq(atom3Ptr->position(),
                                  atom2Ptr->position(), dr3);
         if (needEnergy) {
            dihedralEnergy = interaction().energy(dr1, dr2, dr3, type);
            fraction = (isLocal0 + isLocal1 + isLocal2 + isLocal3)*0.25;
            energy += fraction*dihedralEnergy;
         }
         if (needForce) {
            interaction().force(dr1, dr2, dr2, f1, f2, f3, type);
            if (isLocal0) {
               atom0Ptr->force() += f1;
            }
            if (isLocal1) {
               atom1Ptr->force() -= f1;
               atom1Ptr->force() += f2;
            }
            if (isLocal2) {
               atom2Ptr->force() -= f2;
               atom2Ptr->force() += f3;
            }
            if (isLocal3) {
               atom3Ptr->force() -= f3;
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
   void 
   DihedralPotentialImpl<Interaction>::computeStress(MPI::Intracomm& communicator)
   #else
   void DihedralPotentialImpl<Interaction>::computeStress()
   #endif
   {
      // Do nothing and return if stress is already set.
      if (isStressSet()) return;
 

      Tensor localStress;
      Vector dr1, dr2, dr3;
      Vector f1,  f2, f3;
      GroupIterator<4> iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      Atom*  atom2Ptr;
      Atom*  atom3Ptr;
      int    type;
      int    isLocal0, isLocal1, isLocal2, isLocal3;

      localStress.zero();
      // Iterate over bonds
      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         atom2Ptr = iter->atomPtr(2);
         atom3Ptr = iter->atomPtr(3);
         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         isLocal2 = !(atom2Ptr->isGhost());
         isLocal3 = !(atom3Ptr->isGhost());

         // Calculate stress here

      }

      // Normalize by volume 
      localStress /= boundary().volume();

      // Add localStress from all nodes, set sum on master
      reduceStress(localStress, communicator);
   }

}
#endif
