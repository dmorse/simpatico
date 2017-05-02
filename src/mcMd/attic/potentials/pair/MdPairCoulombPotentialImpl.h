#ifndef SIMP_NOPAIR
#ifndef MCMD_MD_PAIR_COULOMB_POTENTIAL_IMPL_H
#define MCMD_MD_PAIR_COULOMB_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/potentials/pair/McPairPotentialImpl.h>
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
   * Implementation template for an MdPairPotential.
   *
   * \ingroup McMd_Pair_Module
   */
   template <class PairInteraction, class CoulombInteraction>
   class MdPairCoulombPotentialImpl : public MdPairPotential
   {

   public:

      /** 
      * Constructor.
      */
      MdPairCoulombPotentialImpl(System& system);

      /** 
      * Constructor (copied from McPairPotentialImpl)
      */
      MdPairCoulombPotentialImpl(McPairPotentialImpl<PairInteraction,CoulombInteraction>& other);

      /** 
      * Destructor.
      */
      virtual ~MdPairCoulombPotentialImpl();

      /**
      * Read pair potential interaction and pair list blocks.
      * 
      * This method reads the maxBoundary, PairList and pair potential 
      * PairInteraction parameter blocks, in  that order, and initializes an
      * internal PairList. Before calling the PairInteraction::readParameters 
      * method, it passes nAtomType to PairInteraction::setNAtomType().
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
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return repulsive force (< 0 if attractive) over distance
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
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      * \param value  new value of parameter
      */
      void set(std::string name, int i, int j, double value)
      {  pairInteractionPtr_->set(name, i, j, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      */
      double get(std::string name, int i, int j) const
      {  return pairInteractionPtr_->get(name, i, j); }

      /**
      * Return pair interaction class name (e.g., "LJPair").
      */
      virtual std::string interactionClassName() const;

      /**
      * Calculate the total nonBonded pair energy for this System.
      * 
      * Rebuilds the PairList if necessary before calculating energy.
      */
      virtual double energy();

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
      * Compute total nonbonded pressure
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const;

      /**
      * Compute x, y, z nonbonded pressures.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const;

      /**
      * Compute nonbonded stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const;

   protected:

      PairInteraction& pairInteraction() const
      {  return *pairInteractionPtr_; }

      CoulombInteraction& coulombPairInteraction() const
      {  return *coulombPairInteractionPtr_; }

   private:
  
      PairInteraction* pairInteractionPtr_;

      CoulombInteraction* coulombPairInteractionPtr_;

      bool       isCopy_;
 
      template <typename T>
      void computeStressImpl(T& stress) const;

   };

}

#include <mcMd/simulation/System.h> 
#include <mcMd/simulation/Simulation.h> 
#include <mcMd/simulation/stress.h>
#include <mcMd/neighbor/PairIterator.h> 
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
   template <class PairInteraction, class CoulombInteraction>
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::MdPairCoulombPotentialImpl(System& system)
    : MdPairPotential(system),
      pairInteractionPtr_(0),
      coulombPairInteractionPtr_(0),
      isCopy_(false)
   {  
      pairInteractionPtr_ = new PairInteraction; 
      coulombPairInteractionPtr_ = new CoulombInteraction; 
   }

   /* 
   * Constructor, copy from McPairPotentialImpl<PairInteraction,CoulombInteraction>.
   */
   template <class PairInteraction, class CoulombInteraction>
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::MdPairCoulombPotentialImpl(
                         McPairPotentialImpl<PairInteraction,CoulombInteraction>& other)
    : MdPairPotential(other.system()),
      pairInteractionPtr_(&other.pairInteraction()),
      isCopy_(true)
   {}
 
   /* 
   * Destructor. 
   */
   template <class PairInteraction, class CoulombInteraction>
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::~MdPairCoulombPotentialImpl() 
   {
      if (pairInteractionPtr_ && !isCopy_) {
         delete pairInteractionPtr_;
         pairInteractionPtr_ = 0;
      }
      if (coulombInteractionPtr_ && !isCopy_) {
         delete coulombInteractionPtr_;
         coulombInteractionPtr_ = 0;
      }
   }

   template <class PairInteraction, class CoulombInteraction>
   void MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::readParameters(std::istream& in)
   {
      // Read pair potential parameters only if not a copy.
      if (!isCopy_) {
         pairInteraction().setNAtomType(simulation().nAtomType());
         bool nextIndent = false;
         addParamComposite(pairInteraction(), nextIndent);
         pairInteraction().readParameters(in);
      }

      read<Boundary>(in, "maxBoundary", maxBoundary_);
      readParamComposite(in, pairList_);

      // Note: pairlist_ is allocated in MdPairPotential::buildPairList 
      // the first time that function is called. Delaying allocation allows 
      // information about the coulomb rCutoff to be used to choose the 
      // cutoff, because CoulombPotential appears after PairPotential in
      // the parameter file. 
   }

   /*
   * Load internal state from an archive.
   */
   template <class PairInteraction, class CoulombInteraction>
   void 
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::loadParameters(Serializable::IArchive &ar)
   {
      ar >> isCopy_;
      if (!isCopy_) {
         pairInteraction().setNAtomType(simulation().nAtomType());
         bool nextIndent = false;
         addParamComposite(pairInteraction(), nextIndent);
         pairInteraction().loadParameters(ar);
      }
      loadParameter<Boundary>(ar, "maxBoundary", maxBoundary_);
      loadParamComposite(ar, pairList_);
   }

   /*
   * Save internal state to an archive.
   */
   template <class PairInteraction, class CoulombInteraction>
   void MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::save(Serializable::OArchive &ar)
   {
      ar << isCopy_;
      if (!isCopy_) {
         pairInteraction().save(ar);
      }
      ar << maxBoundary_;
      pairList_.save(ar);
   }

   /*
   * Return pair energy for a single pair.
   */
   template <class PairInteraction, class CoulombInteraction> double 
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::energy(double rsq, 
                                        int iAtomType, int jAtomType) const
   { return pairInteraction().energy(rsq, iAtomType, jAtomType); }

   /*
   * Return force / separation for a single pair.
   */
   template <class PairInteraction, class CoulombInteraction> double 
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::forceOverR(double rsq, 
                                        int iAtomType, int jAtomType) const
   {
      if (rsq < pairInteraction().cutoffSq(iAtomType, jAtomType)) { 
         return pairInteraction().forceOverR(rsq, iAtomType, jAtomType); 
      } else {
         return 0.0;
      }
   }

   /*
   * Return maximum cutoff.
   */
   template <class PairInteraction, class CoulombInteraction>
   double MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::maxPairCutoff() const
   {
      double pairCutoff = pairInteraction().maxPairCutoff();
      double coulombCutoff = coulombPairInteraction().rCutoff();
      max = (coulombCutoff > pairCutoff ? coulombCutoff : pairCutoff);
      return max;
   }
 
   /*
   * Return pair interaction class name.
   */
   template <class PairInteraction, class CoulombInteraction>
   std::string MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::interactionClassName() const
   {  return pairInteraction().className(); }

   /* 
   * Return nonBonded Pair interaction energy
   */
   template <class PairInteraction, class CoulombInteraction>
   double MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::energy()
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
      double energy=0;
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position());
         energy += pairInteraction().
                   energy(rsq, atom0Ptr->typeId(), atom1Ptr->typeId());
         energy += coulombPairInteraction().
                   energy(rsq, atom0Ptr->typeId(), atom1Ptr->typeId());
      }

      return energy;
   }
 
   /* 
   * Add nonBonded pair forces to atomic forces.
   */
   template <class PairInteraction, class CoulombInteraction>
   void 
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::addForces()
   {
      // Update PairList if necessary
      if (!isPairListCurrent()) {
         buildPairList();
      }

      PairIterator  iter;
      Vector  force;
      double  rsq, fOverR;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int  type0, type1;

      // Loop over nonbonded neighbor pairs
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position(), force);
         type0 = atom0Ptr->typeId();
         type1 = atom1Ptr->typeId();
         if (rsq < pairInteraction().cutoffSq(type0, type1)) { 
            fOverR = pairInteraction().forceOverR(rsq, type0, type1);
                   + coulombPairInteraction().forceOverR(rsq, type0, type1);
            force *= fOverR;
            atom0Ptr->force() += force;
            atom1Ptr->force() -= force;
         }
      }

   } 

   /* 
   * Add nonBonded pair forces to atomic forces.
   */
   template <class PairInteraction, class CoulombInteraction>
   template <typename T>
   void 
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::computeStressImpl(T& stress) const
   {
      Vector  dr;
      Vector  force;
      double  rsq, fOverR;
      PairIterator iter;
      Atom*  atom1Ptr;
      Atom*  atom0Ptr;
      int  type0, type1;

      setToZero(stress);

      // Loop over nonbonded neighbor pairs
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position(), dr);
         type0 = atom0Ptr->typeId();
         type1 = atom1Ptr->typeId();
         if (rsq < pairInteraction().cutoffSq(type0, type1)) {
            force  = dr;
            fOverR = pairInteraction().forceOverR(rsq, type0, type1);
                   + coulombPairInteraction().forceOverR(rsq, type0, type1);
            force *= forceOverR;
            incrementPairStress(force, dr, stress);
         }
      }

      // Normalize by volume 
      stress /= boundary().volume();
      normalizeStress(stress);
   }

   template <class PairInteraction, class CoulombInteraction>
   void 
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::computeStress(double& stress) const
   {  computeStressImpl(stress); }

   template <class PairInteraction, class CoulombInteraction>
   void 
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::computeStress(Util::Vector& stress) const
   {  computeStressImpl(stress); }

   template <class PairInteraction, class CoulombInteraction>
   void 
   MdPairCoulombPotentialImpl<PairInteraction,CoulombInteraction>::computeStress(Util::Tensor& stress) const
   {  computeStressImpl(stress); }

}
#endif
#endif
