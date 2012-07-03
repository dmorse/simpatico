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
      * block. Before calling Interaction::readParam(), it passes
      * simulation().nDihedralType() to Interaction::setNAtomType().
      */
      virtual void readParam(std::istream& in);

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
      virtual
      double energy(const Vector& R1, const Vector& R2, const Vector& R3,
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
      virtual
      void force(const Vector& R1, const Vector& R2, const Vector& R3,
                 Vector& F1, Vector& F2, Vector& F3, int type) const;

      #if 0
      /**
      * Return pair interaction class name (e.g., "CosineDihedral").
      */
      virtual std::string interactionClassName() const;
      #endif

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
      virtual void addForces();

      /**
      * Add the dihedral forces for all atoms and compute energy.
      */
      virtual void addForces(double& energy);

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
      * Get the total dihedral energy, computed previously by computeEnergy().
      *
      * Call only on master. 
      */
      virtual double energy();

      #if 0
      /**
      * Compute total dihedral pressure.
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const;

      /**
      * Compute x, y, z dihedral pressure components.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const;

      /**
      * Compute dihedral stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const;

      //@}
      #endif

   private:

      /**
      * Total dihedral energy on all processors.
      */
      double energy_;
 
      /**
      * Pointer to Interaction (evaluates the dihedral potential function).
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

#include "DihedralPotential.h"
#include <ddMd/simulation/Simulation.h>
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
   void DihedralPotentialImpl<Interaction>::readParam(std::istream& in)
   {
      readBegin(in,"DihedralPotential");
      bool nextIndent = false;
      readParamComposite(in, interaction(), nextIndent);
      readEnd(in);
   }

   /*
   * Return energy for a single dihedral.
   */
   template <class Interaction>
   inline double DihedralPotentialImpl<Interaction>::
      energy(const Vector& R1, const Vector& R2, const Vector& R3, int typeId) const
   {  return interaction().energy(R1, R2, R3, typeId); }

   /*
   * Return forces for a single dihedral.
   */
   template <class Interaction>
   inline void DihedralPotentialImpl<Interaction>::
      force(const Vector& R1, const Vector& R2, const Vector& R3,
            Vector& F1, Vector& F2, Vector& F3, int typeId) const
   {  interaction().force(R1, R2, R3, F1, F2, F3, typeId); }

   #if 0
   /*
   * Return dihedral potential interaction class name.
   */
   template <class Interaction>
   std::string DihedralPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }
   #endif

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
   void DihedralPotentialImpl<Interaction>::addForces()
   {  addForces(true, false);  }

   /*
   * Increment atomic forces and compute pair energy for this processor.
   */
   template <class Interaction>
   void DihedralPotentialImpl<Interaction>::addForces(double& energy)
   {  energy = addForces(true, true);  }

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
   double DihedralPotentialImpl<Interaction>::energy()
   {  return energy_; } 

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction> double 
   DihedralPotentialImpl<Interaction>::addForces(bool needForce, bool needEnergy)
   {
      // Preconditions
      //if (!storagePtr_->isInitialized()) {
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

      storagePtr_->begin(iter);
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
         boundaryPtr_->distanceSq(atom1Ptr->position(),
                                  atom0Ptr->position(), dr1);
         boundaryPtr_->distanceSq(atom2Ptr->position(),
                                  atom1Ptr->position(), dr2);
         boundaryPtr_->distanceSq(atom3Ptr->position(),
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

   #if 0
   /* 
   * Compute dihedral contribution to stress.
   */
   template <class Interaction>
   template <typename T>
   void DihedralPotentialImpl<Interaction>::computeStressImpl(T& stress) const
   {
      Vector dr1, dr2, f1, f2;
      const Atom *atom0Ptr, *atom1Ptr, *atom2Ptr;

      setToZero(stress);

      // Iterate over all dihedrals in System.
      storagePtr_->begin(iter);
      for ( ; iter.notEnd(); ++iter){
         atom0Ptr = &(iter->atom(0));
         atom1Ptr = &(iter->atom(1));
         atom2Ptr = &(iter->atom(2));

         boundary().distanceSq(atom1Ptr->position(),
                               atom0Ptr->position(), dr1);
         boundary().distanceSq(atom2Ptr->position(),
                               atom1Ptr->position(), dr2);

         // f1 -- along dr1; f2 -- along dr2.
         interaction().force(dr1, dr2,
                           f1, f2, iter->typeId());

         dr1 *= -1.0;
         dr2 *= -1.0;
         incrementPairStress(f1, dr1, stress);
         incrementPairStress(f2, dr2, stress);
      }

      stress /= boundary().volume();
      normalizeStress(stress);
   }

   template <class Interaction>
   void DihedralPotentialImpl<Interaction>::computeStress(double& stress) const
   {  computeStressImpl(stress); }

   template <class Interaction>
   void DihedralPotentialImpl<Interaction>::computeStress(Util::Vector& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Interaction>
   void DihedralPotentialImpl<Interaction>::computeStress(Util::Tensor& stress) 
        const
   {  computeStressImpl(stress); }
   #endif

}
#endif
