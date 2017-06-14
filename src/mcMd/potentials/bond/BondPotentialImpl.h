#ifndef MCMD_BOND_POTENTIAL_IMPL_H
#define MCMD_BOND_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/simulation/SystemInterface.h>

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
   * Implementation template for a BondPotential.
   *
   * \ingroup McMd_Bond_Module
   */
   template <class Interaction>
   class BondPotentialImpl : public BondPotential, private SystemInterface
   {

   public:

      /** 
      * Constructor.
      */
      BondPotentialImpl(System& system);

      /** 
      * Constructor.
      */
      BondPotentialImpl(BondPotentialImpl<Interaction>& other);

      /** 
      * Destructor.
      */
      virtual ~BondPotentialImpl();

      /**
      * Read potential energy.
      * 
      * This method reads the bond potential Interaction parameter
      * block. Before calling Interaction::readParameters(), it passes
      * simulation().nBondType() to Interaction::setNAtomType().
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

      /// \name Interactions Methods
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

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param type   bond type index 
      * \param value  new value of parameter
      */
      void set(std::string name, int type, double value)
      {   interactionPtr_->set(name, type, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param type bond type index 1
      */
      double get(std::string name, int type) const
      {   return interactionPtr_->get(name, type); }

      /**
      * Return pair interaction class name (e.g., "HarmonicBond").
      */
      virtual std::string interactionClassName() const;

      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Calculate the bond energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return bond energy of one atom. 
      */
      double atomEnergy(const Atom& atom) const;

      /**
      * Add the bond forces for all atoms.
      */
      void addForces();

      /**
      * Calculate and store pair energy for this System.
      */
      virtual void computeEnergy();

      /**
      * Compute and store the total nonbonded pressure
      */
      virtual void computeStress();

      //@}

      /**
      * Return bond interaction by reference.
      */
      Interaction& interaction();

      /**
      * Return bond interaction by const reference.
      */
      const Interaction& interaction() const;

   private:
  
      Interaction* interactionPtr_;

      bool isCopy_;

      /**
      * Generic stress calculation, T=Tensor, Vector or double.
      */  
      template <typename T>
      void computeStressImpl(T& stress) const;

   };

}

#include <mcMd/simulation/System.h> 
#include <mcMd/simulation/Simulation.h> 
#include <mcMd/simulation/stress.h>
#include <mcMd/chemistry/getAtomGroups.h>

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
   BondPotentialImpl<Interaction>::BondPotentialImpl(System& system)
    : BondPotential(),
      SystemInterface(system),
      interactionPtr_(0),
      isCopy_(false)
   { interactionPtr_ = new Interaction(); }
 
   /* 
   * Constructor.
   */
   template <class Interaction>
   BondPotentialImpl<Interaction>::BondPotentialImpl(
                         BondPotentialImpl<Interaction>& other)
    : BondPotential(),
      SystemInterface(other.system()),
      interactionPtr_(&other.interaction()),
      isCopy_(true)
   {}
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   BondPotentialImpl<Interaction>::~BondPotentialImpl() 
   {
      if (!isCopy_) {
         delete interactionPtr_;
         interactionPtr_ = 0;
      }
   }

   /* 
   * Read parameters from file.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::readParameters(std::istream &in) 
   {
      // Read only if not a copy.  Do not indent interaction block.
      if (!isCopy_) {
         interaction().setNBondType(simulation().nBondType());
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
   BondPotentialImpl<Interaction>::loadParameters(Serializable::IArchive &ar)
   {
      ar >> isCopy_;
      if (!isCopy_) {
         interaction().setNBondType(simulation().nBondType());
         bool nextIndent = false;
         addParamComposite(interaction(), nextIndent);
         interaction().loadParameters(ar);
      } 
   }

   /*
   * Save internal state to an archive.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::save(Serializable::OArchive &ar)
   {
      ar << isCopy_;
      if (!isCopy_) {
         interaction().save(ar);
      }
   }

   /*
   * Return bond energy for a single pair.
   */
   template <class Interaction>
   double BondPotentialImpl<Interaction>::energy(double rsq, int iBondType) 
      const
   { return interaction().energy(rsq, iBondType); }

   /*
   * Return force / separation for a single bonded pair.
   */
   template <class Interaction>
   double BondPotentialImpl<Interaction>::forceOverR(double rsq, int iBondType)
      const
   { return interaction().forceOverR(rsq, iBondType); }

   /*
   * Return force / separation for a single bonded pair.
   */
   template <class Interaction> double 
   BondPotentialImpl<Interaction>::
      randomBondLength(Random* random, double beta, int bondTypeId) const
   { return interaction().randomBondLength(random, beta, bondTypeId); }

   /*
   * Return bond energy for one Atom. 
   */
   template <class Interaction>
   double BondPotentialImpl<Interaction>::atomEnergy(const Atom &atom) const
   {

      AtomBondArray bonds;
      const Bond* bondPtr;
      double rsq;
      double energy   = 0.0;
      int iBond, nBond;

      getAtomBonds(atom, bonds);
      nBond = bonds.size();
      for (iBond = 0; iBond < nBond; ++iBond) {
         bondPtr = bonds[iBond];
         rsq = boundary().distanceSq(bondPtr->atom(0).position(),
                                     bondPtr->atom(1).position());
         energy += interaction().energy(rsq, bondPtr->typeId());
      }

      return energy;
   }

   /* 
   * Calculate covalent bond energy.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::computeEnergy()
   {
      double energy = 0.0;
      double rsq    = 0.0;
      System::ConstMoleculeIterator molIter;
      Molecule::ConstBondIterator bondIter;

      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (simulation().species(iSpec).nBond() > 0) {
            for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               for (molIter->begin(bondIter); bondIter.notEnd(); ++bondIter) {
                  rsq = boundary().
                        distanceSq( bondIter->atom(0).position(), 
                                    bondIter->atom(1).position());
                  energy += interaction().energy(rsq, bondIter->typeId());
               }
            } 
         }
      }

      energy_.set(energy);
      //return energy;
   }

   /* 
   * Add bonded pair forces to forces array.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::addForces() 
   {
      Vector force;
      double rsq;
      System::MoleculeIterator molIter;
      Molecule::BondIterator bondIter;
      Atom *atom0Ptr, *atom1Ptr;
      int iSpec;

      // Loop over all bonds in system
      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         if (simulation().species(iSpec).nBond() > 0) {
            for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               for (molIter->begin(bondIter); bondIter.notEnd(); ++bondIter) {
                  atom0Ptr = &(bondIter->atom(0));
                  atom1Ptr = &(bondIter->atom(1));
                  rsq = boundary().distanceSq(atom0Ptr->position(), 
                                              atom1Ptr->position(), force);
                  force *= interaction().forceOverR(rsq, bondIter->typeId());
                  atom0Ptr->force() += force;
                  atom1Ptr->force() -= force;
               }
            }
         }
      }
   }

   /* 
   * Generic stress / pressure computation.
   *
   * Allowed types: 
   *    T == Tensor -> compute stress tensor
   *    T == Vector -> compute xx, yy, zz diagonal components
   *    T == double -> compute pressure (P_xx + P_yy + P_zz)/3
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

   /*
   * Compute all short-range pair contributions to stress.
   */
   template <class Interaction>
   void BondPotentialImpl<Interaction>::computeStress()
   {
      Tensor stress;
      computeStressImpl(stress);

      // Set value of Setable<double> energy_ 
      stress_.set(stress);
   }

   template <class Interaction>
   inline Interaction& BondPotentialImpl<Interaction>::interaction()
   { return *interactionPtr_; }

   template <class Interaction>
   inline 
   const Interaction& BondPotentialImpl<Interaction>::interaction() const
   { return *interactionPtr_; }

   /*
   * Return bond potential interaction class name.
   */
   template <class Interaction>
   std::string BondPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }

}
#endif
