#ifndef MCMD_ANGLE_POTENTIAL_IMPL_H
#define MCMD_ANGLE_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/angle/AnglePotential.h>  // base class
#include <mcMd/simulation/SubSystem.h>             // base class
#include <util/global.h>

namespace Util
{
   class Vector;
   class Tensor;
}

namespace McMd
{

   using namespace Util;

   class System;

   /**
   * Implementation template for an AnglePotential.
   *
   * \ingroup McMd_Angle_Module
   */
   template <class Interaction>
   class AnglePotentialImpl : public AnglePotential, public SubSystem
   {

   public:

      /** 
      * Constructor.
      */
      AnglePotentialImpl(System& system);

      /** 
      * Constructor (copied from McAnglePotential)
      */
      AnglePotentialImpl(AnglePotentialImpl<Interaction>& other);

      /** 
      * Destructor.
      */
      virtual ~AnglePotentialImpl();

      /**
      * Read angle potential parameters.
      * 
      * This method reads the angle potential Interaction parameter 
      * block.  Before calling Evalutor::readParameters(), it passes 
      * simulation().nBondType() to Interaction::setNAtomType().
      *
      * \param in input parameter stream.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /// \name Interactions Interface
      //@{

      /**
      * Returns potential energy for one angle.
      *
      * \param cosTheta  cosine of the bend angle.
      * \param type      type of bend angle.
      */
      double energy(double cosTheta, int type) const;
 
      /**
      * Returns forces along two bonds at the angle, for use in MD and stress
      * calculation.
      *
      * \param R1     bond vector from atom 1 to 2.
      * \param R2     bond vector from atom 2 to 3.
      * \param F1     return force along R1 direction.
      * \param F2     return force along R2 direction.
      * \param type   type of angle.
      */
      void force(const Vector& R1, const Vector& R2,
                       Vector& F1, Vector& F2, int type) const;

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name  parameter name
      * \param type  angle type index 
      * \param value new value of parameter
      */
      void set(std::string name, int type, double value)
      {   interactionPtr_->set(name, type, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name  parameter name
      * \param type  angle type index
      * \return parameter value
      */
      double get(std::string name, int type) const
      {   return interactionPtr_->get(name, type); }

      /**
      * Return pair interaction class name (e.g., "CosineAngle").
      */
      virtual std::string interactionClassName() const;

      /**
      * Return angle interaction by reference.
      */
      Interaction& interaction();

      /**
      * Return angle interaction by const reference.
      */
      const Interaction& interaction() const;

      //@}
      /// \name System Energy and Force Calculators
      //@{

      /**
      * Calculate the angle energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return angle energy of one atom. 
      */
      virtual double atomEnergy(const Atom& atom) const;

      /**
      * Calculate the angle energy for one Atom.
      */
      virtual void addForces();

      /**
      * Return total angle potential energy of this System.
      */
      virtual double energy() const;

      /**
      * Compute total angle pressure
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const;

      /**
      * Compute x, y, z angle pressures.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const;

      /**
      * Compute angle stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const;

      //@}

   private:
  
      Interaction* interactionPtr_;

      bool isCopy_;
 
      template <typename T>
      void computeStressImpl(T& stress) const;

   };

}

#include <mcMd/simulation/System.h> 
#include <mcMd/simulation/Simulation.h> 
#include <mcMd/simulation/stress.h>
#include <mcMd/species/Species.h>
#include <util/boundary/Boundary.h> 

#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/accumulators/setToZero.h>

#include <fstream>

namespace McMd
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   template <class Interaction>
   AnglePotentialImpl<Interaction>::AnglePotentialImpl(System& system)
    : AnglePotential(),
      SubSystem(system),
      interactionPtr_(0),
      isCopy_(false)
   {  interactionPtr_ = new Interaction(); }
 
   /* 
   * Constructor, copy from AnglePotentialImpl<Interaction>.
   */
   template <class Interaction>
   AnglePotentialImpl<Interaction>::AnglePotentialImpl(
                         AnglePotentialImpl<Interaction>& other)
    : AnglePotential(),
      SubSystem(other.system()),
      interactionPtr_(&other.interaction()),
      isCopy_(true)
   {}
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   AnglePotentialImpl<Interaction>::~AnglePotentialImpl() 
   {
      if (interactionPtr_ && !isCopy_) {
         delete interactionPtr_; 
         interactionPtr_ = 0;
      }
   }

   /* 
   * Read parameters from file.
   */
   template <class Interaction>
   void AnglePotentialImpl<Interaction>::readParameters(std::istream &in) 
   {
      if (!isCopy_) {
         interaction().setNAngleType(simulation().nAngleType());
         bool nextIndent = false;
         addParamComposite(interaction(), nextIndent);
         interaction().readParameters(in);
      }
   }

   /*
   * Load internal state from an archive.
   */
   template <class Interaction>
   void 
   AnglePotentialImpl<Interaction>::loadParameters(Serializable::IArchive &ar)
   {
      ar >> isCopy_;
      if (!isCopy_) {
         interaction().setNAngleType(simulation().nAngleType());
         bool nextIndent = false;
         addParamComposite(interaction(), nextIndent);
         interaction().loadParameters(ar);
      } 
   }

   /*
   * Save internal state to an archive.
   */
   template <class Interaction>
   void AnglePotentialImpl<Interaction>::save(Serializable::OArchive &ar)
   {
      ar << isCopy_;
      if (!isCopy_) {
         interaction().save(ar);
      }
   }

   /*
   * Return energy for a single angle.
   */
   template <class Interaction>
   double AnglePotentialImpl<Interaction>::energy(double cosTheta, int angleTypeId) 
      const
   {  return interaction().energy(cosTheta, angleTypeId); }

   /*
   * Return forces for a single angle.
   */
   template <class Interaction>
   void AnglePotentialImpl<Interaction>::force(const Vector& R1, const Vector& R2, 
                                          Vector& F1, Vector& F2, int typeId) const
   {  interaction().force(R1, R2, F1, F2, typeId); }

   /*
   * Return angle energy for one Atom. 
   */
   template <class Interaction>
   double AnglePotentialImpl<Interaction>::atomEnergy(const Atom &atom) const
   {
      Species::AtomAngleArray angles;
      const  Angle* anglePtr;
      int    iAngle;
      Vector dr1; // R[1] - R[0]
      Vector dr2; // R[2] - R[1]
      double rsq1, rsq2, cosTheta;
      double energy = 0.0;

      atom.molecule().species().getAtomAngles(atom, angles);
      for (iAngle = 0; iAngle < angles.size(); ++iAngle) {
         anglePtr = angles[iAngle];
         rsq1 = boundary().distanceSq(anglePtr->atom(1).position(),
                               anglePtr->atom(0).position(), dr1);
         rsq2 = boundary().distanceSq(anglePtr->atom(2).position(),
                               anglePtr->atom(1).position(), dr2);
         cosTheta = dr1.dot(dr2) / sqrt(rsq1 * rsq2);
         energy += interaction().energy(cosTheta, anglePtr->typeId());
      }

      return energy;
   }

   /* 
   * Calculate angle energy.
   */
   template <class Interaction>
   double AnglePotentialImpl<Interaction>::energy() const
   {
      Vector dr1; // R[1] - R[0]
      Vector dr2; // R[2] - R[1]
      double rsq1, rsq2, cosTheta;
      double energy = 0.0;
      System::ConstMoleculeIterator  molIter;
      Molecule::ConstAngleIterator angleIter;

      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (simulation().species(iSpec).nAngle() > 0) {
            for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               for (molIter->begin(angleIter); angleIter.notEnd(); ++angleIter){
                  rsq1 = boundary().distanceSq(angleIter->atom(1).position(),
                                        angleIter->atom(0).position(), dr1);
                  rsq2 = boundary().distanceSq(angleIter->atom(2).position(),
                                        angleIter->atom(1).position(), dr2);
                  cosTheta = dr1.dot(dr2) / sqrt(rsq1 * rsq2);
                  energy += interaction().energy(cosTheta, angleIter->typeId());
               }
            }
         }
      }

      return energy;
   }

   /* 
   * Add angle forces to forces array.
   */
   template <class Interaction>
   void AnglePotentialImpl<Interaction>::addForces() 
   {
      Vector dr1, dr2, force1, force2;
      System::MoleculeIterator molIter;
      Molecule::AngleIterator angleIter;
      Atom *atom0Ptr, *atom1Ptr, *atom2Ptr;
      int iSpec;

      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (simulation().species(iSpec).nAngle() > 0) {
            for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               for (molIter->begin(angleIter); angleIter.notEnd(); ++angleIter){
                  atom0Ptr = &(angleIter->atom(0));
                  atom1Ptr = &(angleIter->atom(1));
                  atom2Ptr = &(angleIter->atom(2));
                  boundary().distanceSq(atom1Ptr->position(),
                                        atom0Ptr->position(), dr1);
                  boundary().distanceSq(atom2Ptr->position(),
                                        atom1Ptr->position(), dr2);
                  interaction().force(dr1, dr2, force1, force2,
                                    angleIter->typeId());
                  atom0Ptr->force() += force1;
                  atom1Ptr->force() -= force1;
                  atom1Ptr->force() += force2;
                  atom2Ptr->force() -= force2;
               }
            }
         }
      }
   }

   /* 
   * Compute angle contribution to stress.
   */
   template <class Interaction>
   template <typename T>
   void AnglePotentialImpl<Interaction>::computeStressImpl(T& stress) const
   {
      Vector dr1, dr2, force1, force2;
      const Atom *atom0Ptr, *atom1Ptr, *atom2Ptr;

      setToZero(stress);

      // Iterate over all angles in System.
      System::ConstMoleculeIterator molIter;
      Molecule::ConstAngleIterator  angleIter;
      for (int iSpec = 0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (simulation().species(iSpec).nAngle() > 0) {
            for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               for (molIter->begin(angleIter); angleIter.notEnd(); ++angleIter){
                  atom0Ptr = &(angleIter->atom(0));
                  atom1Ptr = &(angleIter->atom(1));
                  atom2Ptr = &(angleIter->atom(2));

                  boundary().distanceSq(atom1Ptr->position(),
                                        atom0Ptr->position(), dr1);
                  boundary().distanceSq(atom2Ptr->position(),
                                        atom1Ptr->position(), dr2);

                  // force1 -- along dr1; force2 -- along dr2.
                  interaction().force(dr1, dr2,
                                    force1, force2, angleIter->typeId());

                  dr1 *= -1.0;
                  dr2 *= -1.0;
                  incrementPairStress(force1, dr1, stress);
                  incrementPairStress(force2, dr2, stress);
               }
            }
         }
      }

      stress /= boundary().volume();
      normalizeStress(stress);
   }

   template <class Interaction>
   void AnglePotentialImpl<Interaction>::computeStress(double& stress) const
   {  computeStressImpl(stress); }

   template <class Interaction>
   void AnglePotentialImpl<Interaction>::computeStress(Util::Vector& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Interaction>
   void AnglePotentialImpl<Interaction>::computeStress(Util::Tensor& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Interaction>
   inline Interaction& AnglePotentialImpl<Interaction>::interaction()
   { return *interactionPtr_; }

   template <class Interaction>
   inline const Interaction& AnglePotentialImpl<Interaction>::interaction() const
   { return *interactionPtr_; }

   /*
   * Return angle potential interaction class name.
   */
   template <class Interaction>
   std::string AnglePotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }

}
#endif
