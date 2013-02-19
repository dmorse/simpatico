#ifndef DDMD_ANGLE_POTENTIAL_IMPL_H
#define DDMD_ANGLE_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AnglePotential.h" // base class

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
   * Implementation template for a AnglePotential.
   *
   * \ingroup DdMd_Angle_Module
   */
   template <class Interaction>
   class AnglePotentialImpl : public AnglePotential
   {

   public:

      /** 
      * Constructor.
      */
      AnglePotentialImpl(Simulation& simulation);

      /** 
      * Default constructor.
      */
      AnglePotentialImpl();

      /** 
      * Destructor.
      */
      virtual ~AnglePotentialImpl();

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNAngleType(int nAtomType);
  
      /**
      * Read potential energy parameters.
      * 
      * This method reads the angle potential Interaction parameter
      * block. Before calling Interaction::readParameters(), it passes
      * simulation().nAngleType() to Interaction::setNAtomType().
      */
      virtual void readParameters(std::istream& in);

      /// \name Interaction interface
      //@{

      /**
      * Returns potential energy for one angle.
      *
      * \param cosTheta  cosine of the bend angle
      * \param type      type index of angle group
      */
      double angleEnergy(double cosTheta, int type) const;
 
      /**
      * Returns forces along two bonds at the angle, for use in MD and stress
      * calculation.
      *
      * \param R1     bond vector from atom 1 to 2
      * \param R2     bond vector from atom 2 to 3
      * \param F1     return force along R1 direction
      * \param F2     return force along R2 direction
      * \param type   type index for angle group
      */
      void angleForce(const Vector& R1, const Vector& R2,
                      Vector& F1, Vector& F2, int type) const;

      /**
      * Modify an angle parameter, identified by a string.
      *
      * \param name  parameter variable name
      * \param type  type index for angle group
      * \param value  new value of parameter
      */
      void set(std::string name, int type, double value)
      {  interactionPtr_->set(name, type, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name  parameter variable name
      * \param type  type index of angle group
      */
      double get(std::string name, int type) const
      {  return interactionPtr_->get(name, type); }

      /**
      * Return pair interaction class name (e.g., "CosineAngle").
      */
      virtual std::string interactionClassName() const;

      /**
      * Return angle interaction by const reference.
      */
      const Interaction& interaction() const;

      /**
      * Return angle interaction by reference.
      */
      Interaction& interaction();

      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Compute angle forces for atoms.
      *
      * Adds angle forces to existing Atom forces.
      */
      virtual void computeForces();

      /**
      * Compute the total angle energy for all processors
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
      * Pointer to Interaction (evaluates the angle potential function).
      */ 
      Interaction* interactionPtr_;

   };

}

#include "AnglePotential.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/GroupStorage.h>
#include <ddMd/storage/GroupIterator.h>

#include <util/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/global.h>

#include <fstream>

namespace DdMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   template <class Interaction>
   AnglePotentialImpl<Interaction>::AnglePotentialImpl(Simulation& simulation)
    : AnglePotential(simulation),
      interactionPtr_(0)
   {  interactionPtr_ = new Interaction(); }
 
   /* 
   * Default constructor.
   */
   template <class Interaction>
   AnglePotentialImpl<Interaction>::AnglePotentialImpl()
    : AnglePotential(),
      interactionPtr_(0)
   {  interactionPtr_ = new Interaction(); }
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   AnglePotentialImpl<Interaction>::~AnglePotentialImpl() 
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
   void AnglePotentialImpl<Interaction>::setNAngleType(int nAngleType)
   {  interaction().setNAngleType(nAngleType); }

   /*
   * Read angle interaction parameters.
   */
   template <class Interaction>
   void AnglePotentialImpl<Interaction>::readParameters(std::istream& in)
   {
      bool nextIndent = false;
      addParamComposite(interaction(), nextIndent);
      interaction().readParameters(in);
   }

   /*
   * Return energy for a single angle.
   */
   template <class Interaction>
   double 
   AnglePotentialImpl<Interaction>::angleEnergy(double cosTheta, int angleTypeId) 
      const
   {  return interaction().energy(cosTheta, angleTypeId); }

   /*
   * Return forces for a single angle.
   */
   template <class Interaction> void 
   AnglePotentialImpl<Interaction>::angleForce(const Vector& R1, 
                                               const Vector& R2, 
                                               Vector& F1, Vector& F2, 
                                               int typeId) const
   {  interaction().force(R1, R2, F1, F2, typeId); }

   /*
   * Return angle potential interaction class name.
   */
   template <class Interaction>
   std::string AnglePotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }

   /**
   * Get Interaction by reference.
   */
   template <class Interaction>
   inline Interaction& AnglePotentialImpl<Interaction>::interaction()
   {  return *interactionPtr_; }

   /**
   * Get Interaction by const reference.
   */
   template <class Interaction>
   inline const Interaction& AnglePotentialImpl<Interaction>::interaction() const
   {  return *interactionPtr_; }

   /*
   * Increment atomic forces, without calculating energy.
   */
   template <class Interaction>
   void AnglePotentialImpl<Interaction>::computeForces()
   {  
      Vector dr1; // R[1] - R[0]
      Vector dr2; // R[2] - R[1]
      Vector f1;  // d(energy)/d(dr1)
      Vector f2;  // d(energy)/d(dr2)
      GroupIterator<3> iter;
      Atom* atom0Ptr;
      Atom* atom1Ptr;
      Atom* atom2Ptr;
      int   type;

      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         atom2Ptr = iter->atomPtr(2);
         // Calculate minimimum image separations
         boundary().distanceSq(atom1Ptr->position(),
                               atom0Ptr->position(), dr1);
         boundary().distanceSq(atom2Ptr->position(),
                               atom1Ptr->position(), dr2);
         interaction().force(dr1, dr2, f1, f2, type);
         if (!atom0Ptr->isGhost()) {
            atom0Ptr->force() += f1;
         }
         if (!atom1Ptr->isGhost()) {
            atom1Ptr->force() -= f1;
            atom1Ptr->force() += f2;
         }
         if (!atom2Ptr->isGhost()) {
            atom2Ptr->force() -= f2;
         }
      }
   }

   /*
   * Compute total angle energy on all processors.
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void 
   AnglePotentialImpl<Interaction>::computeEnergy(MPI::Intracomm& communicator)
   #else
   void AnglePotentialImpl<Interaction>::computeEnergy()
   #endif
   {

      // Do nothing and return if energy is already set.
      if (isEnergySet()) return;
 
      Vector dr1; // R[1] - R[0]
      Vector dr2; // R[2] - R[1]
      double rsq1, rsq2, cosTheta;
      double angleEnergy;
      double localEnergy = 0.0;
      double fraction;
      double third = 1.0/3.0;
      GroupIterator<3> iter;
      Atom* atom0Ptr;
      Atom* atom1Ptr;
      Atom* atom2Ptr;
      int   type, isLocal0, isLocal1, isLocal2;

      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         type = iter->typeId();
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         atom2Ptr = iter->atomPtr(2);
         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         isLocal2 = !(atom2Ptr->isGhost());
         // Calculate minimimum image separations
         rsq1 = boundary().distanceSq(atom1Ptr->position(),
                                      atom0Ptr->position(), dr1);
         rsq2 = boundary().distanceSq(atom2Ptr->position(),
                                      atom1Ptr->position(), dr2);
         cosTheta = dr1.dot(dr2) / sqrt(rsq1 * rsq2);
         angleEnergy = interaction().energy(cosTheta, type);
         fraction = (isLocal0 + isLocal1 + isLocal2)*third;
         localEnergy += fraction*angleEnergy;
      }

      // Add localEnergy from all nodes, set sum on master.
      reduceEnergy(localEnergy, communicator);
   }

   /*
   * Compute total pair stress (Call on all processors).
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void AnglePotentialImpl<Interaction>::computeStress(MPI::Intracomm& communicator)
   #else
   void AnglePotentialImpl<Interaction>::computeStress()
   #endif
   {
      // Do nothing and return if stress is already set.
      if (isStressSet()) return;
 
      Tensor localStress;
      Vector dr1, dr2;
      Vector f1, f2;
      double factor;
      double prefactor = -1.0/3.0;
      GroupIterator<3> iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      Atom*  atom2Ptr;
      int    type;
      int    isLocal0, isLocal1, isLocal2;

      // Iterate over angle groups
      localStress.zero();
      storage().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         atom2Ptr = iter->atomPtr(2);
         type = iter->typeId();
         boundary().distanceSq(atom1Ptr->position(),
                                  atom0Ptr->position(), dr1);
         boundary().distanceSq(atom2Ptr->position(),
                                  atom1Ptr->position(), dr2);

         // Calculate derivatives f1, f2 of energy with respect to dr1, dr2
         interaction().force(dr1, dr2, f1, f2, type);

         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         isLocal2 = !(atom2Ptr->isGhost());
         factor = prefactor*(isLocal0 + isLocal1 + isLocal2);

         // Increment localStress tensor
         dr1 *= factor;
         dr2 *= factor;
         incrementPairStress(f1, dr1, localStress);
         incrementPairStress(f2, dr2, localStress);
      }

      // Normalize by volume 
      localStress /= boundary().volume();

      // Add localEnergy from all nodes, set sum on master.
      reduceStress(localStress, communicator);
   }

}
#endif
