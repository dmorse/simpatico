#ifndef MCMD_NOPAIR
#ifndef MC_PAIR_INTERACTION_IMPL_H
#define MC_PAIR_INTERACTION_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/boundary/Boundary.h>
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
   * \ingroup Pair_Module
   */
   template <class Evaluator>
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
      * Reads maxBoundary and pair potential Evaluator blocks.
      * 
      * This method reads the maxBoundary and Evaluator parameter blocks,
      * and then allocates memory for an internal CellList. It passes
      * simulation().nAtomType() as an argument Evaluator::setNAtomType() 
      * before calling Evaluator::readParam().
      */
      virtual void readParam(std::istream& in);

      /// \name Energy, Force, Stress Evaluators
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
      * Return pair evaluator class name (e.g., "LJPair").
      */
      virtual std::string evaluatorClassName() const;

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
      * Return total nonbonded pair potential energy of this System.
      */
      double energy();

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
      * Compute stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const;

      //@}

      /**
      * Return reference to underlying pair evaluator.
      */
      Evaluator& evaluator();

      /**
      * Return reference to underlying pair evaluator.
      */
      const Evaluator& evaluator() const;

   private:
  
      Evaluator evaluator_;
 
      template <typename T>
      void computeStressImpl(T& stress) const;

   };

}

#include <mcMd/simulation/System.h> 
#include <mcMd/simulation/Simulation.h> 
#include <mcMd/simulation/stress.h>
#include <mcMd/neighbor/PairIterator.h> 
#include <mcMd/boundary/Boundary.h> 

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
   McPairPotentialImpl<Evaluator>::McPairPotentialImpl(System& system)
    : McPairPotential(system)
   {}
 
   /* 
   * Destructor. 
   */
   template <class Evaluator>
   McPairPotentialImpl<Evaluator>::~McPairPotentialImpl() 
   {}

   /* 
   * Read parameters from file.
   */
   template <class Evaluator>
   void McPairPotentialImpl<Evaluator>::readParam(std::istream &in) 
   {
      readBegin(in, "McPairPotential");
     
      // Must setNAtomTypes in evaluator before calling readParam.
      evaluator().setNAtomType(simulation().nAtomType());

      // Read potential energy parameters with no indent or brackets.
      bool nextIndent = false;
      readParamComposite(in, evaluator(), nextIndent);

      // Read maxBoundary (needed to allocate memory for cell list).
      read<Boundary>(in, "maxBoundary", maxBoundary_);

      // Allocate memory for the CellList.
      cellList_.allocate(simulation().atomCapacity(), maxBoundary_, 
                         maxPairCutoff());

      readEnd(in);
   }
  
   /*
   * Return pair energy for a single pair.
   */
   template <class Evaluator>
   double McPairPotentialImpl<Evaluator>::energy(double rsq, int iAtomType, int jAtomType) const
   { return evaluator().energy(rsq, iAtomType, jAtomType); }

   /*
   * Return force / separation for a single pair.
   */
   template <class Evaluator>
   double McPairPotentialImpl<Evaluator>::forceOverR(double rsq, int iAtomType, int jAtomType) const
   { return evaluator().forceOverR(rsq, iAtomType, jAtomType); }

   /*
   * Return maximum cutoff.
   */
   template <class Evaluator>
   double McPairPotentialImpl<Evaluator>::maxPairCutoff() const
   { return evaluator().maxPairCutoff(); }

   /* 
   * Return nonbonded pair energy for one Atom.
   */
   template <class Evaluator>
   double McPairPotentialImpl<Evaluator>::atomEnergy(const Atom &atom) const
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
               energy += evaluator().
                         energy(rsq, atom.typeId(), jAtomPtr->typeId());
            }
         }
      } 
      return energy;
   }

   /* 
   * Return nonbonded pair potential energy for one Molecule.
   */
   template <class Evaluator>
   double McPairPotentialImpl<Evaluator>::moleculeEnergy(const Molecule &molecule) const
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
                     energy += evaluator().energy(rsq, iAtomPtr->typeId(), 
                                                           jAtomPtr->typeId());
                  }
               }
            }
         }
      } 
      return energy;
   }

   /*
   * Return total nonbonded pair potential energy for System.
   */
   template <class Evaluator>
   double McPairPotentialImpl<Evaluator>::energy()
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
                     energy += evaluator().energy(rsq, iAtomPtr->typeId(), 
                                                  jAtomPtr->typeId());
                  }

               }

            } // secondary atoms

         } // primary atoms

      } // cells

      return energy;
   }

   /*
   * Return nonbonded pair stress or pressure (implementation template).
   */
   template <class Evaluator>
   template <typename T>
   void McPairPotentialImpl<Evaluator>::computeStressImpl(T& stress) const
   {
      Vector force, dr;
      double rsq;
      const Atom *atom0Ptr, *atom1Ptr;
      int nNeighbor, nInCell, ia, ja;

      setToZero(stress);

      // Loop over cells
      for (int ic=0; ic < cellList_.totCells(); ++ic) {

         // Get array of neighbors
         cellList_.getCellNeighbors(ic, neighbors_, nInCell);
         nNeighbor = neighbors_.size();
  
         // Loop over primary atoms in this cell.
         for (ia = 0; ia < nInCell; ++ia) {
            atom0Ptr = neighbors_[ia];
          
            // Loop over secondary atoms in the neighboring cells.
            for (ja = 0; ja < nNeighbor; ++ja) {
               atom1Ptr = neighbors_[ja];
     
               // Count each pair only once.
               if (atom1Ptr->id() > atom0Ptr->id()) {

                  // Exclude masked pairs.
                  if (!atom0Ptr->mask().isMasked(*atom1Ptr)) {

                     // Evaluate force.
                     rsq = boundary().distanceSq(atom0Ptr->position(),
                                                 atom1Ptr->position(), dr);
                     force = dr;
                     force *= evaluator().forceOverR(rsq,
                              atom0Ptr->typeId(), atom1Ptr->typeId());

                     incrementPairStress(force, dr, stress);

                  }

               }

            } // secondary atoms

         } // primary atoms

      } // cells

      // Normalize by volume.
      stress /= boundary().volume();
      normalizeStress(stress);

   }


   template <class Evaluator>
   void McPairPotentialImpl<Evaluator>::computeStress(double& stress) const
   {  computeStressImpl(stress); }


   template <class Evaluator>
   void McPairPotentialImpl<Evaluator>::computeStress(Util::Vector& stress) 
        const
   {  computeStressImpl(stress); }


   template <class Evaluator>
   void McPairPotentialImpl<Evaluator>::computeStress(Util::Tensor& stress) 
        const
   {  computeStressImpl(stress); }


   template <class Evaluator>
   inline Evaluator& McPairPotentialImpl<Evaluator>::evaluator()
   { return evaluator_; }

   template <class Evaluator>
   inline const Evaluator& McPairPotentialImpl<Evaluator>::evaluator() const
   { return evaluator_; }

   /*
   * Return pair evaluator class name.
   */
   template <class Evaluator>
   std::string McPairPotentialImpl<Evaluator>::evaluatorClassName() const
   {  return evaluator().className(); }

}
#endif
#endif
