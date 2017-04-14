#ifndef MCMD_MD_PAIR_POTENTIAL_IMPL_H
#define MCMD_MD_PAIR_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/potentials/pair/McPairPotentialImpl.h>
#include <mcMd/simulation/stress.h>
#include <util/space/Tensor.h>
#include <util/misc/Setable.h>
#include <util/global.h>

namespace Util
{
   class Vector;
}

namespace McMd
{

   using namespace Util;

   class System;

   /**
   * Implementation template for an MdPairPotential.
   *
   * \ingroup McMd_Pair_Module
   */
   template <class Interaction>
   class MdPairPotentialImpl : public MdPairPotential
   {

   public:

      /**
      * Constructor.
      */
      MdPairPotentialImpl(System& system);

      /**
      * Constructor (copied from McPairPotentialImpl)
      */
      MdPairPotentialImpl(McPairPotentialImpl<Interaction>& other);

      /**
      * Destructor.
      */
      virtual ~MdPairPotentialImpl();

      /**
      * Read pair potential interaction and pair list blocks.
      *
      * This method reads the pair potential Interaction parameter 
      * and PairList blocks, and initializes an internal PairList. 
      * Before calling the Interaction::readParameters method, it 
      * passes nAtomType to Interaction::setNAtomType().
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

      /// \name Pair Interaction Interface
      //@{

      /**
      * Return pair energy for a single pair.
      *
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return energy of pair
      */
      virtual
      double energy(double rsq, int iAtomType, int jAtomType) const;

      /**
      * Return force / separation for a single pair.
      *
      * \param rsq  square distance between atoms in pair
      * \param iAtomType  atom type index of 1st atom
      * \param jAtomType  atom type index of 2nd atom
      * \return  repulsive force (< 0 if attractive) over distance
      */
      virtual
      double forceOverR(double rsq, int iAtomType, int jAtomType) const;

      /**
      * Return maximum cutoff.
      */
      virtual double maxPairCutoff() const;

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name  parameter name
      * \param i  type index of first atom
      * \param j  type index of first atom
      * \param value  new value of parameter
      */
      void set(std::string name, int i, int j, double value)
      {  interactionPtr_->set(name, i, j, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name  parameter name
      * \param i  type index of first atom
      * \param j  type index of first atom
      */
      double get(std::string name, int i, int j) const
      {  return interactionPtr_->get(name, i, j); }

      /**
      * Return pair interaction class name (e.g., "LJPair").
      */
      virtual std::string interactionClassName() const;

      //@}
      /// \name Global force and energy calculators
      //@{

      /**
      * Calculate non-bonded pair forces for all atoms in this System.
      *
      * Adds non-bonded pair forces to the current values of the
      * forces for all atoms in this system. Before calculating
      * forces, the method checks if the pair list is current,
      * and rebuilds it if necessary.
      */
      virtual void addForces();

      /**
      * Calculate and store pair energy for this System.
      *
      * Rebuilds the PairList if necessary before calculating energy.
      */
      virtual void computeEnergy();

      /**
      * Compute and store the total nonbonded pressure
      */
      virtual void computeStress();

      //@}

   protected:

      Interaction& interaction() const
      {  return *interactionPtr_; }

      /*
      * Generalized stress computation.
      */
      template <typename T>
      void computeStressImpl(T& stress);

   private:

      Interaction* interactionPtr_;

      bool       isCopy_;

   };

}

#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/stress.h>
#include <mcMd/neighbor/PairIterator.h>
#include <util/boundary/Boundary.h>

#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/accumulators/setToZero.h>

#include <fstream>

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   template <class Interaction>
   MdPairPotentialImpl<Interaction>::MdPairPotentialImpl(System& system)
    : MdPairPotential(system),
      interactionPtr_(0),
      isCopy_(false)
   {  interactionPtr_ = new Interaction; }

   /*
   * Constructor, copy from McPairPotentialImpl<Interaction>.
   */
   template <class Interaction>
   MdPairPotentialImpl<Interaction>::MdPairPotentialImpl(
                         McPairPotentialImpl<Interaction>& other)
    : MdPairPotential(other.system()),
      interactionPtr_(&other.interaction()),
      isCopy_(true)
   {}

   /*
   * Destructor.
   */
   template <class Interaction>
   MdPairPotentialImpl<Interaction>::~MdPairPotentialImpl()
   {
      if (interactionPtr_ && !isCopy_) {
         delete interactionPtr_;
         interactionPtr_ = 0;
      }
   }

   template <class Interaction>
   void MdPairPotentialImpl<Interaction>::readParameters(std::istream& in)
   {
      // Read pair potential parameters only if not a copy.
      if (!isCopy_) {
         interaction().setNAtomType(simulation().nAtomType());
         bool nextIndent = false;
         addParamComposite(interaction(), nextIndent);
         interaction().readParameters(in);
      }

      // Initialize the PairList
      readParamComposite(in, pairList_);
      double cutoff = interaction().maxPairCutoff();
      pairList_.initialize(simulation().atomCapacity(), cutoff);
   }

   /*
   * Load internal state from an archive.
   */
   template <class Interaction>
   void
   MdPairPotentialImpl<Interaction>::loadParameters(Serializable::IArchive &ar)
   {
      ar >> isCopy_;
      if (!isCopy_) {
         interaction().setNAtomType(simulation().nAtomType());
         bool nextIndent = false;
         addParamComposite(interaction(), nextIndent);
         interaction().loadParameters(ar);
      }
      loadParamComposite(ar, pairList_);
   }

   /*
   * Save internal state to an archive.
   */
   template <class Interaction>
   void MdPairPotentialImpl<Interaction>::save(Serializable::OArchive &ar)
   {
      ar << isCopy_;
      if (!isCopy_) {
         interaction().save(ar);
      }
      pairList_.save(ar);
   }

   /*
   * Return pair energy for a single pair.
   */
   template <class Interaction> double
   MdPairPotentialImpl<Interaction>::energy(double rsq,
                                        int iAtomType, int jAtomType) const
   {  return interaction().energy(rsq, iAtomType, jAtomType); }

   /*
   * Return force / separation for a single pair.
   */
   template <class Interaction> double
   MdPairPotentialImpl<Interaction>::forceOverR(double rsq,
                                        int iAtomType, int jAtomType) const
   {
      if (rsq < interaction().cutoffSq(iAtomType, jAtomType)) {
         return interaction().forceOverR(rsq, iAtomType, jAtomType);
      } else {
         return 0.0;
      }
   }

   /*
   * Return maximum cutoff.
   */
   template <class Interaction>
   double MdPairPotentialImpl<Interaction>::maxPairCutoff() const
   { return interaction().maxPairCutoff(); }

   /*
   * Return pair interaction class name.
   */
   template <class Interaction>
   std::string MdPairPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }

   /*
   * Add nonBonded pair forces to atomic forces.
   */
   template <class Interaction>
   void MdPairPotentialImpl<Interaction>::addForces()
   {
      // Update PairList if necessary
      if (!isPairListCurrent()) {
         buildPairList();
      }

      PairIterator iter;
      Vector       force;
      double       rsq;
      Atom        *atom0Ptr;
      Atom        *atom1Ptr;
      int          type0, type1;

      // Loop over nonbonded neighbor pairs
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position(),
                          force);
         type0 = atom0Ptr->typeId();
         type1 = atom1Ptr->typeId();
         if (rsq < interaction().cutoffSq(type0, type1)) {
            force *= interaction().forceOverR(rsq, type0, type1);
            atom0Ptr->force() += force;
            atom1Ptr->force() -= force;
         }
      }

   }

   /*
   * Compute and store all short-range pair energy components.
   */
   template <class Interaction>
   void MdPairPotentialImpl<Interaction>::computeEnergy()
   {
      // Update PairList if necessary
      if (!isPairListCurrent()) {
         buildPairList();
      }

      // Loop over pairs
      PairIterator iter;
      Atom *atom0Ptr;
      Atom *atom1Ptr;
      double rsq;
      double energy = 0.0;
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position());
         energy += interaction().
                   energy(rsq, atom0Ptr->typeId(), atom1Ptr->typeId());
      }

      // Set value of Setable<double> energy_ 
      energy_.set(energy);
   }

   /*
   * Compute all short-range pair contributions to stress.
   */
   template <class Interaction>
   template <typename T>
   void MdPairPotentialImpl<Interaction>::computeStressImpl(T& stress)
   {
      // Update PairList if necessary
      if (!isPairListCurrent()) {
         buildPairList();
      }

      Vector dr;
      Vector force;
      double rsq;
      PairIterator iter;
      Atom* atom1Ptr;
      Atom* atom0Ptr;
      int type0, type1;

      // Set all elements of stress tensor to zero.
      setToZero(stress);

      // Loop over nonbonded neighbor pairs
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position(), dr);
         type0 = atom0Ptr->typeId();
         type1 = atom1Ptr->typeId();
         if (rsq < interaction().cutoffSq(type0, type1)) {
            force  = dr;
            force *= interaction().forceOverR(rsq, type0, type1);
            incrementPairStress(force, dr, stress);
         }
      }

      // Normalize by volume
      stress /= boundary().volume();
      normalizeStress(stress);
   }

   /*
   * Compute all short-range pair contributions to stress.
   */
   template <class Interaction>
   void MdPairPotentialImpl<Interaction>::computeStress()
   {
      Tensor stress;
      computeStressImpl(stress);

      // Set value of Setable<double> stress_
      stress_.set(stress);
   }

}
#endif
