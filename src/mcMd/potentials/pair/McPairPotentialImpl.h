#ifndef MCMD_MC_PAIR_INTERACTION_IMPL_H
#define MCMD_MC_PAIR_INTERACTION_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/pair/McPairPotential.h>
#include <util/boundary/Boundary.h>
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
   * Implementation template for an McPairPotential.
   *
   * \ingroup McMd_Pair_Module
   */
   template <class Interaction>
   class McPairPotentialImpl : public McPairPotential
   {

   public:

      /** 
      * Constructor.
      */
      McPairPotentialImpl(System& system);

      /** 
      * Destructor.
      */
      virtual ~McPairPotentialImpl();

      /**
      * Reads the pair potential Interaction blocks.
      * 
      * This method reads the Interaction parameter blocks, and initializes
      * an internal CellList. It passes the value simulation().nAtomType() 
      * to the Interaction::setNAtomType() function before calling 
      * Interaction::readParameters().
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

      /// \name Interaction interface
      //@{

      /**
      * Return pair energy for a single pair.
      */
      virtual double energy(double rsq, int iAtomType, int jAtomType) const;

      /**
      * Return force / separation for a single pair.
      */
      virtual double forceOverR(double rsq, int iAtomType, int jAtomType) const;

      /**
      * Return maximum cutoff.
      */
      virtual double maxPairCutoff() const;

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      * \param value  parameter value
      */
      void set(std::string name, int i, int j, double value)
      {  interaction_.set(name, i, j, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      */
      double get(std::string name, int i, int j) const
      {  return interaction_.get(name, i, j); }

      /**
      * Return pair interaction class name (e.g., "LJPair").
      */
      virtual std::string interactionClassName() const;

      /**
      * Return reference to underlying pair interaction.
      */
      Interaction& interaction();

      /**
      * Return reference to underlying pair interaction.
      */
      const Interaction& interaction() const;

      //@}
      /// \name Energy and stress calculators
      //@{

      /**
      * Calculate the nonbonded pair energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return nonbonded pair potential energy of atom
      */
      double atomEnergy(const Atom& atom) const;

      /**
      * Calculate the nonbonded pair energy for an entire Molecule.
      *
      * \param  molecule Molecule object of interest
      * \return nonbonded pair potential energy of atom
      */
      double moleculeEnergy(const Molecule& molecule) const;

      /**
      * Compute and store nonbonded pair energy of this System.
      *
      * If energy is already known (set), this function does nothing.
      */
      virtual void computeEnergy();

      /**
      * Compute and store total short-range pair stress.
      *
      * If stress is already known (set), this function does nothing.
      */
      virtual void computeStress();

      //@}

   private:
 
      /**
      * Pair interaction object (e.g., Interaction == LJPair)
      */ 
      Interaction interaction_;

      /**
      * Generalized stress computation, for variable type T.
      *
      * Allowed types: T =, double, Vector or Tensor 
      * For T == Tensor, computes stress tensor
      * For T == Vector, computes xx, yy, zz components
      * For T == double, computes pressure
      */
      template <typename T>
      void computeStressImpl(T& stress);

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
   template <class Interaction>
   McPairPotentialImpl<Interaction>::McPairPotentialImpl(System& system)
    : McPairPotential(system)
   {}
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   McPairPotentialImpl<Interaction>::~McPairPotentialImpl() 
   {}

   /* 
   * Read parameters from file.
   */
   template <class Interaction>
   void McPairPotentialImpl<Interaction>::readParameters(std::istream &in) 
   {
      interaction().setNAtomType(simulation().nAtomType());

      // Read potential energy parameters with no indent or brackets.
      bool nextIndent = false;
      addParamComposite(interaction(), nextIndent);
      interaction().readParameters(in);

      // Set atom capacity and allocate memory in the CellList.
      cellList_.setAtomCapacity(simulation().atomCapacity());
   }

   /*
   * Load internal state from an archive.
   */
   template <class Interaction>
   void 
   McPairPotentialImpl<Interaction>::loadParameters(Serializable::IArchive &ar)
   {
      interaction().setNAtomType(simulation().nAtomType());

      // Read potential energy parameters with no indent or brackets.
      bool nextIndent = false;
      addParamComposite(interaction(), nextIndent);
      interaction().loadParameters(ar);

      // Allocate memory for the CellList.
      cellList_.setAtomCapacity(simulation().atomCapacity());
   }

   /*
   * Save internal state to an archive.
   */
   template <class Interaction>
   void McPairPotentialImpl<Interaction>::save(Serializable::OArchive &ar)
   {
      interaction().save(ar);
   }

   /*
   * Return pair energy for a single pair.
   */
   template <class Interaction>
   double 
   McPairPotentialImpl<Interaction>::energy(double rsq, int iAtomType, int jAtomType) 
   const
   {  return interaction().energy(rsq, iAtomType, jAtomType); }

   /*
   * Return force / separation for a single pair.
   */
   template <class Interaction>
   double 
   McPairPotentialImpl<Interaction>::forceOverR(double rsq, int iAtomType, int jAtomType) 
   const
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
   double McPairPotentialImpl<Interaction>::maxPairCutoff() const
   {  return interaction().maxPairCutoff(); }

   /* 
   * Return nonbonded pair energy for one Atom.
   */
   template <class Interaction>
   double McPairPotentialImpl<Interaction>::atomEnergy(const Atom &atom) const
   {
      Atom   *jAtomPtr;
      double  energy;
      double  rsq;
      int     j, jId, nNeighbor;
      int     id = atom.id();

      // Get array of neighbors
      cellList_.getNeighbors(atom.position(), neighbors_);
      nNeighbor = neighbors_.size();

      // Loop over neighboring atoms
      energy = 0.0;
      for (j = 0; j < nNeighbor; ++j) {
         jAtomPtr = neighbors_[j];
         jId      = jAtomPtr->id();

         // Check if atoms are the same
         if (jId != id) {

            // Check if atoms are bonded
            if (!atom.mask().isMasked(*jAtomPtr)) {
               rsq = boundary().
                     distanceSq(atom.position(), jAtomPtr->position());
               energy += interaction().
                         energy(rsq, atom.typeId(), jAtomPtr->typeId());
            }
         }
      } 
      return energy;
   }

   /* 
   * Return nonbonded pair potential energy for one Molecule.
   */
   template <class Interaction>
   double McPairPotentialImpl<Interaction>::moleculeEnergy(const Molecule &molecule) const
   {
      const Atom* iAtomPtr;
      const Atom* jAtomPtr;
      double  energy;
      double  rsq;
      int     i, iId, j, jId, nNeighbor;

      energy = 0.0;
      for (i = 0; i < molecule.nAtom(); ++i) { 
         iAtomPtr = &molecule.atom(i);
         iId      = iAtomPtr->id();

         // Get array of neighbors
         cellList_.getNeighbors(iAtomPtr->position(), neighbors_);
         nNeighbor = neighbors_.size();
   
         // Loop over neighboring atoms
         for (j = 0; j < nNeighbor; ++j) {
            jAtomPtr = neighbors_[j];
            jId      = jAtomPtr->id();
   
            // Check if atoms are the same
            if (jId != iId) {
   
               // Check if atoms are bonded
               if (!(iAtomPtr->mask().isMasked(*jAtomPtr))) {

                  if ( (&iAtomPtr->molecule() != &jAtomPtr->molecule()) ||
                       (iId < jId) ) {
                     rsq = boundary().distanceSq(iAtomPtr->position(), 
                                                 jAtomPtr->position());
                     energy += interaction().energy(rsq, iAtomPtr->typeId(), 
                                                           jAtomPtr->typeId());
                  }
               }
            }
         }
      } 
      return energy;
   }

   /*
   * Compute and store total nonbonded pair potential energy for System.
   */
   template <class Interaction>
   void McPairPotentialImpl<Interaction>::computeEnergy()
   {
      Atom  *iAtomPtr, *jAtomPtr;
      double energy;
      double  rsq;
      int    nNeighbor, nInCell;
      int    i, j, iId, jId, ic;

      // Loop over cells
      energy = 0.0;
      for (ic=0; ic < cellList_.totCells(); ++ic) {

         // Get array of neighbors
         cellList_.getCellNeighbors(ic, neighbors_, nInCell);
         nNeighbor = neighbors_.size();
  
         // Loop over primary atoms in this cell
         for (i = 0; i < nInCell; ++i) {
            iAtomPtr = neighbors_[i];
            iId      = iAtomPtr->id();
          
            // Loop over secondary atoms in this and neighboring cells
            for (j = 0; j < nNeighbor; ++j) {
               jAtomPtr = neighbors_[j];
               jId      = jAtomPtr->id();
     
               // Count each pair only once 
               if (jId > iId) {

                  // Exclude masked pairs
                  if (!iAtomPtr->mask().isMasked(*jAtomPtr)) {

                     rsq = boundary().distanceSq(iAtomPtr->position(), 
                                                 jAtomPtr->position());
                     energy += interaction().energy(rsq, iAtomPtr->typeId(), 
                                                  jAtomPtr->typeId());
                  }

               }

            } // secondary atoms

         } // primary atoms

      } // cells

      // Store local variable energy in Setable<double> class member energy_
      energy_.set(energy);
   }

   /*
   * Compute nonbonded pair stress.
   */
   template <class Interaction>
   template <typename T>
   void McPairPotentialImpl<Interaction>::computeStressImpl(T& stress)
   {
      Vector dr;
      Vector force;
      double rsq;
      const Atom *atom0Ptr, *atom1Ptr;
      int nNeighbor, nInCell, ia, ja, type0, type1;

      setToZero(stress);

      // Loop over cells
      for (int ic=0; ic < cellList_.totCells(); ++ic) {

         // Get array of neighbors
         cellList_.getCellNeighbors(ic, neighbors_, nInCell);
         nNeighbor = neighbors_.size();
  
         // Loop over primary atoms in this cell.
         for (ia = 0; ia < nInCell; ++ia) {
            atom0Ptr = neighbors_[ia];
            type0 = atom0Ptr->typeId();
          
            // Loop over secondary atoms in the neighboring cells.
            for (ja = 0; ja < nNeighbor; ++ja) {
               atom1Ptr = neighbors_[ja];
               type1 = atom1Ptr->typeId();
     
               // Count each pair only once.
               if (atom1Ptr->id() > atom0Ptr->id()) {

                  // Exclude masked pairs.
                  if (!atom0Ptr->mask().isMasked(*atom1Ptr)) {
                     rsq = boundary().distanceSq(atom0Ptr->position(),
                                                 atom1Ptr->position(), dr);
                     if (rsq < interaction().cutoffSq(type0, type1)) {
                        force = dr;
                        force *= interaction().forceOverR(rsq, type0, type1);
                        incrementPairStress(force, dr, stress);
                     }

                  }

               }

            } // secondary atoms

         } // primary atoms

      } // cells

      // Normalize by volume.
      stress /= boundary().volume();
      normalizeStress(stress);
   }

   /*
   * Compute nonbonded pair stress.
   */
   template <class Interaction>
   void McPairPotentialImpl<Interaction>::computeStress()
   {
      Tensor stress;
      computeStressImpl(stress);

      stress_.set(stress);
   }

   /*
   * Get non-const reference to Interaction.
   */
   template <class Interaction>
   inline Interaction& McPairPotentialImpl<Interaction>::interaction()
   {  return interaction_; }

   /*
   * Get const reference to Interaction.
   */
   template <class Interaction>
   inline const Interaction& McPairPotentialImpl<Interaction>::interaction() const
   {  return interaction_; }

   /*
   * Return pair interaction class name.
   */
   template <class Interaction>
   std::string McPairPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }

}
#endif
