#ifdef  INTER_DIHEDRAL
#ifndef DIHEDRAL_POTENTIAL_IMPL_H
#define DIHEDRAL_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/dihedral/DihedralPotential.h>  // base class
#include <mcMd/simulation/SubSystem.h>                   // base class
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
   * Implementation template for an DihedralPotential.
   *
   * \ingroup Dihedral_Module
   */
   template <class Evaluator>
   class DihedralPotentialImpl : public DihedralPotential, public SubSystem
   {

   public:

      /** 
      * Constructor.
      */
      DihedralPotentialImpl(System& system);

      /** 
      * Constructor (copied from DihedralPotential)
      */
      DihedralPotentialImpl(DihedralPotentialImpl<Evaluator>& other);

      /** 
      * Destructor.
      */
      virtual ~DihedralPotentialImpl();

      /**
      * Read dihedral potential parameters.
      * 
      * This method reads the dihedral potential Evaluator parameter 
      * block.  Before calling Evalutor::readParam(), it passes 
      * simulation().nBondType() to Evaluator::setNAtomType().
      *
      * \param in input parameter stream.
      */
      virtual void readParam(std::istream& in);

      /// \name Evaluator Interface.
      //@{

      /**
      * Returns potential energy for one dihedral.
      *
      *     1   3    4
      *     o   o----o
      *      \ /
      *       o 
      *       2 
      *
      * \param R1     bond vector from atom 1 to 2.
      * \param R2     bond vector from atom 2 to 3.
      * \param R3     bond vector from atom 3 to 4.
      * \param type   type of dihedral.
      */
      virtual
      double energy(const Vector& R1, const Vector& R2, const Vector& R3,
                    int type) const;
 
      /**
      * Returns derivatives of energy with respect to bond vectors forming the
      * dihedral group.
      *
      * \param R1     bond vector from atom 1 to 2.
      * \param R2     bond vector from atom 2 to 3.
      * \param R3     bond vector from atom 3 to 4.
      * \param F1     return force along R1 direction.
      * \param F2     return force along R2 direction.
      * \param F3     return force along R2 direction.
      * \param type   type of dihedral.
      */
      virtual
      void force(const Vector& R1, const Vector& R2, const Vector& R3,
                 Vector& F1, Vector& F2, Vector& F3, int type) const;

      /**
      * Return pair evaluator class name (e.g., "CosineDihedral").
      */
      virtual std::string evaluatorClassName() const;

      //@}
      /// \name Global Energy, Force, Stress methods.
      //@{

      /**
      * Compute and return the dihedral energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return dihedral energy of one atom.
      */
      virtual double atomEnergy(const Atom& atom) const;

      /**
      * Compute and return total dihedral potential energy of this System.
      */
      virtual double energy() const;

      /**
      * Add dihedral forces to all atomic forces.
      */
      virtual void addForces();

      /**
      * Compute total dihedral pressure
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const;

      /**
      * Compute x, y, z dihedral pressures.
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

      /**
      * Return dihedral potential evaluator by reference.
      */
      Evaluator& evaluator();

      /**
      * Return dihedral potential evaluator by const reference.
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
   DihedralPotentialImpl<Evaluator>::DihedralPotentialImpl(System& system)
    : DihedralPotential(),
      SubSystem(system),
      evaluatorPtr_(0),
      isCopy_(false)
   {  evaluatorPtr_ = new Evaluator(); }
 
   /* 
   * Constructor, copy from DihedralPotentialImpl<Evaluator>.
   */
   template <class Evaluator>
   DihedralPotentialImpl<Evaluator>::DihedralPotentialImpl(
                         DihedralPotentialImpl<Evaluator>& other)
    : DihedralPotential(),
      SubSystem(other.system()),
      evaluatorPtr_(&other.evaluator()),
      isCopy_(true)
   {}
 
   /* 
   * Destructor. 
   */
   template <class Evaluator>
   DihedralPotentialImpl<Evaluator>::~DihedralPotentialImpl() 
   {}

   /* 
   * Read parameters from file.
   */
   template <class Evaluator>
   void DihedralPotentialImpl<Evaluator>::readParam(std::istream &in) 
   {
      if (!isCopy_) {
         readBegin(in, "DihedralPotential");
         evaluator().setNDihedralType(simulation().nDihedralType());
         bool nextIndent = false;
         readParamComposite(in, evaluator(), nextIndent);
         readEnd(in);
      }
   }

   /*
   * Return energy for a single dihedral.
   */
   template <class Evaluator>
   inline double DihedralPotentialImpl<Evaluator>::
      energy(const Vector& R1, const Vector& R2, const Vector& R3, int typeId) const
   {  return evaluator().energy(R1, R2, R3, typeId); }

   /*
   * Return forces for a single dihedral.
   */
   template <class Evaluator>
   inline void DihedralPotentialImpl<Evaluator>::
      force(const Vector& R1, const Vector& R2, const Vector& R3,
            Vector& F1, Vector& F2, Vector& F3, int typeId) const
   {  evaluator().force(R1, R2, R3, F1, F2, F3, typeId); }

   /*
   * Return dihedral energy for one Atom.
   */
   template <class Evaluator>
   double DihedralPotentialImpl<Evaluator>::atomEnergy(const Atom &atom) 
   const
   {
      Species::AtomDihedralArray dihedrals;
      const  Dihedral* dihedralPtr;
      int    iDihedral;
      Vector dr1; // R[1] - R[0]
      Vector dr2; // R[2] - R[1]
      Vector dr3; // R[3] - R[2]
      double energy = 0.0;

      atom.molecule().species().getAtomDihedrals(atom, dihedrals);
      for (iDihedral = 0; iDihedral < dihedrals.size(); ++iDihedral) {
         dihedralPtr = dihedrals[iDihedral];
         boundary().distanceSq(dihedralPtr->atom(1).position(),
                               dihedralPtr->atom(0).position(), dr1);
         boundary().distanceSq(dihedralPtr->atom(2).position(),
                               dihedralPtr->atom(1).position(), dr2);
         boundary().distanceSq(dihedralPtr->atom(3).position(),
                               dihedralPtr->atom(2).position(), dr3);
         energy += evaluator().energy(dr1, dr2, dr3,
                                             dihedralPtr->typeId());
      }
      return energy;
   }

   /*
   * Compute total dihedral energy for a System.
   */
   template <class Evaluator>
   double DihedralPotentialImpl<Evaluator>::energy() const
   {
      Vector dr1; // R[1] - R[0]
      Vector dr2; // R[2] - R[1]
      Vector dr3; // R[3] - R[2]
      double energy = 0.0;
      System::ConstMoleculeIterator  molIter;
      Molecule::ConstDihedralIterator dihedralIter;

      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (simulation().species(iSpec).nDihedral() > 0) {
            for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               molIter->begin(dihedralIter); 
               for ( ; dihedralIter.notEnd(); ++dihedralIter){
                  boundary().distanceSq(dihedralIter->atom(1).position(),
                                        dihedralIter->atom(0).position(), dr1);
                  boundary().distanceSq(dihedralIter->atom(2).position(),
                                        dihedralIter->atom(1).position(), dr2);
                  boundary().distanceSq(dihedralIter->atom(3).position(),
                                        dihedralIter->atom(2).position(), dr3);
                  energy += evaluator().
                            energy(dr1, dr2, dr3, dihedralIter->typeId());
               }
            }
         }
      }

      return energy;
   }

   /* 
   * Add dihedral forces to atomic forces.
   */
   template <class Evaluator>
   void DihedralPotentialImpl<Evaluator>::addForces() 
   {
      Vector dr1, dr2, dr3, force1, force2, force3;
      Atom *atom0Ptr, *atom1Ptr, *atom2Ptr, *atom3Ptr;
      int iSpec;

      System::MoleculeIterator molIter;
      Molecule::DihedralIterator dihedralIter;
      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (simulation().species(iSpec).nDihedral() > 0) {
            for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               molIter->begin(dihedralIter); 
               for ( ; dihedralIter.notEnd(); ++dihedralIter) {

                  atom0Ptr = &(dihedralIter->atom(0));
                  atom1Ptr = &(dihedralIter->atom(1));
                  atom2Ptr = &(dihedralIter->atom(2));
                  atom3Ptr = &(dihedralIter->atom(3));

                  boundary().distanceSq(atom1Ptr->position(),
                                        atom0Ptr->position(), dr1);
                  boundary().distanceSq(atom2Ptr->position(),
                                        atom1Ptr->position(), dr2);
                  boundary().distanceSq(atom3Ptr->position(),
                                        atom2Ptr->position(), dr3);
                  evaluator().force(dr1, dr2, dr3,
                     force1, force2, force3, dihedralIter->typeId());

                  atom0Ptr->force() += force1;
                  atom1Ptr->force() -= force1;
                  atom1Ptr->force() += force2;
                  atom2Ptr->force() -= force2;
                  atom2Ptr->force() += force3;
                  atom3Ptr->force() -= force3;
               }
            }
         }
      }
   }

   /* 
   * Add dihedral contribution to stress.
   */
   template <class Evaluator>
   template <typename T>
   void DihedralPotentialImpl<Evaluator>::computeStressImpl(T& stress) const
   {
      Vector dr1, dr2, dr3, force1, force2, force3;

      setToZero(stress);

      // Iterate over all dihedrals in System.
      System::ConstMoleculeIterator  molIter;
      Molecule::ConstDihedralIterator dihedralIter;
      for (int iSpec = 0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (simulation().species(iSpec).nDihedral() > 0) {
            for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               molIter->begin(dihedralIter);
               for ( ; dihedralIter.notEnd(); ++dihedralIter){

                  boundary().distanceSq(dihedralIter->atom(1).position(),
                                        dihedralIter->atom(0).position(), dr1);
                  boundary().distanceSq(dihedralIter->atom(2).position(),
                                        dihedralIter->atom(1).position(), dr2);
                  boundary().distanceSq(dihedralIter->atom(3).position(),
                                        dihedralIter->atom(2).position(), dr3);

                  evaluator().force(dr1, dr2, dr3, force1, force2, force3, 
                                    dihedralIter->typeId());

                  dr1 *= -1.0;
                  dr2 *= -1.0;
                  dr3 *= -1.0;
                  incrementPairStress(force1, dr1, stress);
                  incrementPairStress(force2, dr2, stress);
                  incrementPairStress(force3, dr3, stress);
               }
            }
         }
      }

      stress /= boundary().volume();
      normalizeStress(stress);
   }
   template <class Evaluator>
   void DihedralPotentialImpl<Evaluator>::computeStress(double& stress) const
   {  computeStressImpl(stress); }

   template <class Evaluator>
   void DihedralPotentialImpl<Evaluator>::computeStress(Util::Vector& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Evaluator>
   void DihedralPotentialImpl<Evaluator>::computeStress(Util::Tensor& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Evaluator>
   inline Evaluator& DihedralPotentialImpl<Evaluator>::evaluator()
   { return *evaluatorPtr_; }

   template <class Evaluator>
   inline const Evaluator& DihedralPotentialImpl<Evaluator>::evaluator() const
   { return *evaluatorPtr_; }

   /*
   * Return dihedral potential evaluator class name.
   */
   template <class Evaluator>
   std::string DihedralPotentialImpl<Evaluator>::evaluatorClassName() const
   {  return evaluator().className(); }

}
#endif
#endif
