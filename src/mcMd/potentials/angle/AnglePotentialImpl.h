#ifdef  INTER_ANGLE
#ifndef ANGLE_POTENTIAL_IMPL_H
#define ANGLE_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
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
   * \ingroup Angle_Module
   */
   template <class Evaluator>
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
      AnglePotentialImpl(AnglePotentialImpl<Evaluator>& other);

      /** 
      * Destructor.
      */
      virtual ~AnglePotentialImpl();

      /**
      * Read angle potential parameters.
      * 
      * This method reads the angle potential Evaluator parameter 
      * block.  Before calling Evalutor::readParam(), it passes 
      * simulation().nBondType() to Evaluator::setNAtomType().
      *
      * \param in input parameter stream.
      */
      virtual void readParam(std::istream& in);

      /// \name Energy, Force, Stress Evaluators
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
      * Return pair evaluator class name (e.g., "CosineAngle").
      */
      virtual std::string evaluatorClassName() const;

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

      /**
      * Return angle evaluator by reference.
      */
      Evaluator& evaluator();

      /**
      * Return angle evaluator by const reference.
      */
      const Evaluator& evaluator() const;

   private:
  
      Evaluator* evaluatorPtr_;

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
   template <class Evaluator>
   AnglePotentialImpl<Evaluator>::AnglePotentialImpl(System& system)
    : AnglePotential(),
      SubSystem(system),
      evaluatorPtr_(0),
      isCopy_(false)
   {  evaluatorPtr_ = new Evaluator(); }
 
   /* 
   * Constructor, copy from AnglePotentialImpl<Evaluator>.
   */
   template <class Evaluator>
   AnglePotentialImpl<Evaluator>::AnglePotentialImpl(
                         AnglePotentialImpl<Evaluator>& other)
    : AnglePotential(),
      SubSystem(other.system()),
      evaluatorPtr_(&other.evaluator()),
      isCopy_(true)
   {}
 
   /* 
   * Destructor. 
   */
   template <class Evaluator>
   AnglePotentialImpl<Evaluator>::~AnglePotentialImpl() 
   {
      if (evaluatorPtr_ && !isCopy_) {
         delete evaluatorPtr_; 
         evaluatorPtr_ = 0;
      }
   }

   /* 
   * Read parameters from file.
   */
   template <class Evaluator>
   void AnglePotentialImpl<Evaluator>::readParam(std::istream &in) 
   {
      if (!isCopy_) {
         readBegin(in, "AnglePotential");
         evaluator().setNAngleType(simulation().nAngleType());
         bool nextIndent = false;
         readParamComposite(in, evaluator(), nextIndent);
         readEnd(in);
      }
   }

   /*
   * Return energy for a single angle.
   */
   template <class Evaluator>
   double AnglePotentialImpl<Evaluator>::energy(double cosTheta, int angleTypeId) 
      const
   {  return evaluator().energy(cosTheta, angleTypeId); }

   /*
   * Return forces for a single angle.
   */
   template <class Evaluator>
   void AnglePotentialImpl<Evaluator>::force(const Vector& R1, const Vector& R2, 
                                          Vector& F1, Vector& F2, int typeId) const
   {  evaluator().force(R1, R2, F1, F2, typeId); }

   /*
   * Return angle energy for one Atom. 
   */
   template <class Evaluator>
   double AnglePotentialImpl<Evaluator>::atomEnergy(const Atom &atom) const
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
         energy += evaluator().energy(cosTheta, anglePtr->typeId());
      }

      return energy;
   }

   /* 
   * Calculate angle energy.
   */
   template <class Evaluator>
   double AnglePotentialImpl<Evaluator>::energy() const
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
                  energy += evaluator().energy(cosTheta, angleIter->typeId());
               }
            }
         }
      }

      return energy;
   }

   /* 
   * Add angle forces to forces array.
   */
   template <class Evaluator>
   void AnglePotentialImpl<Evaluator>::addForces() 
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
                  evaluator().force(dr1, dr2, force1, force2,
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
   template <class Evaluator>
   template <typename T>
   void AnglePotentialImpl<Evaluator>::computeStressImpl(T& stress) const
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
                  evaluator().force(dr1, dr2,
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

   template <class Evaluator>
   void AnglePotentialImpl<Evaluator>::computeStress(double& stress) const
   {  computeStressImpl(stress); }

   template <class Evaluator>
   void AnglePotentialImpl<Evaluator>::computeStress(Util::Vector& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Evaluator>
   void AnglePotentialImpl<Evaluator>::computeStress(Util::Tensor& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Evaluator>
   inline Evaluator& AnglePotentialImpl<Evaluator>::evaluator()
   { return *evaluatorPtr_; }

   template <class Evaluator>
   inline const Evaluator& AnglePotentialImpl<Evaluator>::evaluator() const
   { return *evaluatorPtr_; }

   /*
   * Return angle potential evaluator class name.
   */
   template <class Evaluator>
   std::string AnglePotentialImpl<Evaluator>::evaluatorClassName() const
   {  return evaluator().className(); }

}
#endif
#endif
