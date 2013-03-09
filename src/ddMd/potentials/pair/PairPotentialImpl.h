#ifndef DDMD_PAIR_POTENTIAL_IMPL_H
#define DDMD_PAIR_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/potentials/pair/PairPotential.h>
#include <util/space/Tensor.h>
#include <util/global.h>

//namespace Util {
//   class Vector;
//}

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Implementation template for a PairPotential.
   *
   * \ingroup DdMd_Pair_Module
   */
   template <class Interaction>
   class PairPotentialImpl : public PairPotential
   {

   public:

      /** 
      * Constructor.
      */
      PairPotentialImpl(Simulation& simulation);

      /** 
      * Default constructor (for unit testing).
      */
      PairPotentialImpl();

      /** 
      * Destructor.
      */
      virtual ~PairPotentialImpl();

      /**
      * Read pair potential interaction and pair list blocks.
      * 
      * This method reads the maxBoundary, PairList and pair potential 
      * Interaction parameter blocks, in  that order, and initializes an
      * internal PairList. Before calling the Interaction::readParam method,
      * it passes simulation().nAtomType() to Interaction::setNAtomType().
      *
      * \param in input parameter stream.
      */
      virtual void readParameters(std::istream& in);

      /// \name Interaction interface
      //@{

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNAtomType(int nAtomType);
  
      /**
      * Return energy for a single pair.
      * 
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return energy of pair
      */
      virtual 
      double pairEnergy(double rsq, int iAtomType, int jAtomType) const;

      /**
      * Return force / separation for a single pair.
      *
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return repulsive force (< 0 if attractive) over distance
      */
      virtual 
      double pairForceOverR(double rsq, int iAtomType, int jAtomType) const;

      /**
      * Modify a pair interaction parameter, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      * \param value  new value of parameter
      */
      void set(std::string name, int i, int j, double value)
      {  interactionPtr_->set(name, i, j, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      */
      double get(std::string name, int i, int j) const
      {  return interactionPtr_->get(name, i, j); }

      /**
      * Return maximum cutoff.
      */
      virtual double maxPairCutoff() const;

      /**
      * Return pair interaction class name (e.g., "LJPair").
      */
      virtual std::string interactionClassName() const;

      /**
      * Return underlying pair interaction object by const reference.
      */
      const Interaction& interaction() const
      {  return *interactionPtr_; }

      /**
      * Return underlying pair interaction object by reference.
      */
      Interaction& interaction()
      {  return *interactionPtr_; }

      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Calculate non-bonded pair forces for all atoms in this Simulation.
      *
      * Adds non-bonded pair forces to the current values of the
      * forces for all atoms in this simulation. Before calculating 
      * forces, the method checks if the pair list is current, and
      * rebuilds it if necessary.
      */
      virtual void computeForces();

      /**
      * Compute the total nonBonded pair energy for all processors
      * 
      * Call on all processors.
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator);
      #else
      virtual void computeEnergy();
      #endif

      /**
      * Compute the total nonBonded stress for all processors
      * 
      * Call on all processors.
      */
      #ifdef UTIL_MPI
      virtual void computeStress(MPI::Intracomm& communicator);
      #else
      virtual void computeStress();
      #endif

      /**
      * Compute total pair energies for all processors
      * Compute nonbonded forces and sress for all processors
      * 
      * Call on all processors.
      */
      #ifdef UTIL_MPI
      virtual void computePairEnergies(MPI::Intracomm& communicator);
      #else
      virtual void computePairEnergies();
      #endif

      /**
      * Compute nonbonded forces and sress for all processors
      * 
      * Call on all processors.
      */
      #ifdef UTIL_MPI
      virtual void computeForcesAndStress(MPI::Intracomm& communicator);
      #else
      virtual void computeForcesAndStress();
      #endif

      //@}

   private:

      /**
      * Pointer to pair interaction object.
      */ 
      Interaction* interactionPtr_;

      // Number of atom types.
      int nAtomType_;

      /**
      * Calculate atomic pair energy, using PairList.
      */
      double energyList();

      /**
      * Calculate atomic pair forces, using PairList.
      */
      void computeForcesList();

      /**
      * Calculate atomic pair forces and/or pair potential energy.
      */
      double energyCell();

      /**
      * Calculate atomic pair forces and/or pair potential energy.
      */
      void computeForcesCell();

      /**
      * Calculate atomic pair energy, using N^2 loop.
      * 
      * Use an O(N^2) double loop over all atoms.
      */
      double energyNSq();

      /**
      * Calculate atomic pair forces, using N^2 loop.
      */
      void computeForcesNSq();

   };

}

#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/neighbor/PairIterator.h>
#include <ddMd/communicate/Domain.h>

#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/accumulators/setToZero.h>
#include <util/global.h>

#include <fstream>

namespace DdMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   template <class Interaction>
   PairPotentialImpl<Interaction>::PairPotentialImpl(Simulation& simulation)
    : PairPotential(simulation),
      interactionPtr_(0)
   {  
      interactionPtr_ = new Interaction;
      nAtomType_ =  simulation.nAtomType();
   }
 
   /* 
   * Default constructor.
   */
   template <class Interaction>
   PairPotentialImpl<Interaction>::PairPotentialImpl()
    : PairPotential(),
      interactionPtr_(0)
   {  interactionPtr_ = new Interaction; }
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   PairPotentialImpl<Interaction>::~PairPotentialImpl() 
   {
      if (interactionPtr_) {
         delete interactionPtr_;
         interactionPtr_ = 0;
      }
   }

   template <class Interaction>
   void PairPotentialImpl<Interaction>::readParameters(std::istream& in)
   {
      bool nextIndent = false;
      addParamComposite(interaction(), nextIndent);
      interaction().readParameters(in);
      readPairListParam(in);
   }

   /**
   * Set the maximum number of atom types.
   */
   template <class Interaction>
   void PairPotentialImpl<Interaction>::setNAtomType(int nAtomType)
   {  interaction().setNAtomType(nAtomType); }

   /*
   * Return pair energy for a single pair.
   */
   template <class Interaction> double 
   PairPotentialImpl<Interaction>::pairEnergy(double rsq, 
                                        int iAtomType, int jAtomType) const
   { return interaction().energy(rsq, iAtomType, jAtomType); }

   /*
   * Return force / separation for a single pair.
   */
   template <class Interaction> double 
   PairPotentialImpl<Interaction>::pairForceOverR(double rsq, 
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
   double PairPotentialImpl<Interaction>::maxPairCutoff() const
   { return interaction().maxPairCutoff(); }

   /*
   * Return pair interaction class name.
   */
   template <class Interaction>
   std::string PairPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }

   /*
   * Increment atomic forces, without calculating energy.
   */
   template <class Interaction>
   void PairPotentialImpl<Interaction>::computeForces()
   {  
       if (methodId() == 0) {
          computeForcesList(); 
       } else
       if (methodId() == 1) {
          computeForcesCell(); 
       } else {
          computeForcesNSq(); 
       }
   }

   /*
   * Compute total pair energy on all processors.
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void 
   PairPotentialImpl<Interaction>::computeEnergy(MPI::Intracomm& communicator)
   #else
   void PairPotentialImpl<Interaction>::computeEnergy()
   #endif
   {
      // Do nothing (return) if energy is already set.
      if (isEnergySet()) return;
 
      double localEnergy = 0.0; 
      if (methodId() == 0) {
         localEnergy = energyList(); 
      } else 
      if (methodId() == 1) {
         localEnergy = energyCell(); 
      } else {
         localEnergy = energyNSq(); 
      }

      // Add localEnergy from all nodes, set energy to sum on master.
      reduceEnergy(localEnergy, communicator);
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction>
   double PairPotentialImpl<Interaction>::energyList()
   {
      Vector f;
      double rsq;
      double energy = 0.0;
      PairIterator iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int    type0, type1;
      if (reverseUpdateFlag()) {
         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            assert(!atom0Ptr->isGhost());
            type0 = atom0Ptr->typeId();
            type1 = atom1Ptr->typeId();
            f.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = f.square();
            energy += interactionPtr_->energy(rsq, type0, type1);
         }
      } else {
         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            assert(!atom0Ptr->isGhost());
            type0 = atom0Ptr->typeId();
            type1 = atom1Ptr->typeId();
            f.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = f.square();
            if (!atom1Ptr->isGhost()) {
               energy += interactionPtr_->energy(rsq, type0, type1);
            } else {
               energy += 0.5*interactionPtr_->energy(rsq, type0, type1);
            }
         }
      }
      return energy;
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction>
   void PairPotentialImpl<Interaction>::computeForcesList()
   {
      Vector f;
      double rsq;
      PairIterator iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int    type0, type1;

      if (reverseUpdateFlag()) {

         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            f.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = f.square();
            type0 = atom0Ptr->typeId();
            type1 = atom1Ptr->typeId();
            if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
               f *= interactionPtr_->forceOverR(rsq, type0, type1);
               atom0Ptr->force() += f;
               atom1Ptr->force() -= f;
            }
         }

      } else {

         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            f.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = f.square();
            type0 = atom0Ptr->typeId();
            type1 = atom1Ptr->typeId();
            if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
               f *= interactionPtr_->forceOverR(rsq, type0, type1);
               atom0Ptr->force() += f;
               if (!atom1Ptr->isGhost()) {
                  atom1Ptr->force() -= f;
               }
            }
         }

      }
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction>
   double PairPotentialImpl<Interaction>::energyCell()
   {
      // Find all neighbors (cell list)
      Cell::NeighborArray neighbors;
      Vector f;
      double rsq;
      double energy = 0.0;
      Atom*  atomPtr0;
      Atom*  atomPtr1;
      const Cell*  cellPtr;
      int    type0, type1, na, nn, i, j;

      // Iterate over linked list of local cells.
      cellPtr = cellList_.begin();
      while (cellPtr) {
         cellPtr->getNeighbors(neighbors, reverseUpdateFlag());
         na = cellPtr->nAtom();
         nn = neighbors.size();
         for (i = 0; i < na; ++i) {
            atomPtr0 = neighbors[i];
            type0 = atomPtr0->typeId();

            // Loop over atoms in this cell
            for (j = 0; j < na; ++j) {
               atomPtr1 = neighbors[j];
               type1 = atomPtr1->typeId();
               if (atomPtr1 > atomPtr0) {
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  energy += interactionPtr_->energy(rsq, type0, type1);
               }
            }

            // Loop over atoms in neighboring cells.
            if (reverseUpdateFlag()) {
               for (j = na; j < nn; ++j) {
                  atomPtr1 = neighbors[j];
                  type1 = atomPtr1->typeId();
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  energy += interactionPtr_->energy(rsq, type0, type1);
               }
            } else {
               for (j = na; j < nn; ++j) {
                  atomPtr1 = neighbors[j];
                  type1 = atomPtr1->typeId();
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  if (!atomPtr1->isGhost()) {
                     energy += interactionPtr_->energy(rsq, type0, type1);
                  } else {
                     energy += 0.5*interactionPtr_->energy(rsq, type0, type1);
                  }
               }
            }

         }
         cellPtr = cellPtr->nextCellPtr();
      } // while (cellPtr) 
      return energy;
   }

   /*
   * Increment atomic forces using Cell List (private).
   */
   template <class Interaction>
   void PairPotentialImpl<Interaction>::computeForcesCell()
   {
      // Find all neighbors (cell list)
      Cell::NeighborArray neighbors;
      Vector f;
      double rsq;
      Atom*  atomPtr0;
      Atom*  atomPtr1;
      const Cell* cellPtr;
      int    type0, type1, na, nn, i, j;

      // Iterate over local cells.
      cellPtr = cellList_.begin();
      while (cellPtr) {
         cellPtr->getNeighbors(neighbors);
         na = cellPtr->nAtom();
         nn = neighbors.size();
         for (i = 0; i < na; ++i) {
            atomPtr0 = neighbors[i];
            type0 = atomPtr0->typeId();
            // Loop over atoms in this cell
            for (j = 0; j < na; ++j) {
               atomPtr1 = neighbors[j];
               type1 = atomPtr1->typeId();
               if (atomPtr1 > atomPtr0) {
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
                     f *= interactionPtr_->forceOverR(rsq, type0, type1);
                     atomPtr0->force() += f;
                     atomPtr1->force() -= f;
                  }
               }
            }

            // Loop over atoms in neighboring cells.
            if (reverseUpdateFlag()) {
               for (j = na; j < nn; ++j) {
                  atomPtr1 = neighbors[j];
                  type1 = atomPtr1->typeId();
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
                     f *= interactionPtr_->forceOverR(rsq, type0, type1);
                     atomPtr0->force() += f;
                     atomPtr1->force() -= f;
                  }
               }
            } else {
               for (j = na; j < nn; ++j) {
                  atomPtr1 = neighbors[j];
                  type1 = atomPtr1->typeId();
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
                     f *= interactionPtr_->forceOverR(rsq, type0, type1);
                     atomPtr0->force() += f;
                     if (!atomPtr1->isGhost()) {
                        atomPtr1->force() -= f;
                     }
                  }
               }
            }

         }
         cellPtr = cellPtr->nextCellPtr();
      } // while (cellPtr) 
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction>
   double PairPotentialImpl<Interaction>::energyNSq()
   {
      Vector f;
      double rsq;
      double energy = 0.0;
      AtomIterator  atomIter0, atomIter1;
      GhostIterator ghostIter;
      int id0, id1, type0, type1;

      // Iterate over local atom 0
      storage().begin(atomIter0);
      for ( ; atomIter0.notEnd(); ++atomIter0) {
         id0 = atomIter0->id();
         type0 = atomIter0->typeId();

         // Iterate over local atom 1
         storage().begin(atomIter1);
         for ( ; atomIter1.notEnd(); ++atomIter1) {
            id1 = atomIter1->id();
            if (id0 < id1) {
               if (!atomIter0->mask().isMasked(id1)) {
                  f.subtract(atomIter0->position(), atomIter1->position());
                  rsq = f.square();
                  type1 = atomIter1->typeId();
                  energy += interactionPtr_->energy(rsq, type0, type1);
               }
            }
         }

         // Iterate over ghost atoms
         storage().begin(ghostIter);
         if (reverseUpdateFlag()) {
            for ( ; ghostIter.notEnd(); ++ghostIter) {
               id1 = ghostIter->id();
               if (id0 < id1) {
                  if (!atomIter0->mask().isMasked(id1)) {
                     f.subtract(atomIter0->position(), ghostIter->position());
                     rsq = f.square();
                     type1 = ghostIter->typeId();
                     energy += interactionPtr_->energy(rsq, type0, type1);
                  }
               }
            }
         } else {
            for ( ; ghostIter.notEnd(); ++ghostIter) {
               id1 = ghostIter->id();
               if (!atomIter0->mask().isMasked(id1)) {
                  f.subtract(atomIter0->position(), ghostIter->position());
                  rsq = f.square();
                  type1 = ghostIter->typeId();
                  energy += 0.5*interactionPtr_->energy(rsq, type0, type1);
               }
            }
         }

      }
      return energy;
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction>
   void PairPotentialImpl<Interaction>::computeForcesNSq()
   {
      Vector f;
      double rsq;
      AtomIterator  atomIter0, atomIter1;
      GhostIterator ghostIter;
      int           type0, type1, id0, id1;

      // Iterate over atom 0
      storage().begin(atomIter0);
      for ( ; atomIter0.notEnd(); ++atomIter0) {
         id0 = atomIter0->id();
         type0 = atomIter0->typeId();

         // Iterate over local atom 1
         storage().begin(atomIter1);
         for ( ; atomIter1.notEnd(); ++atomIter1) {
            id1 = atomIter1->id();
            if (id0 < id1) {
               if (!atomIter0->mask().isMasked(id1)) {
                  // Set f = r0 - r1, separation between atoms
                  f.subtract(atomIter0->position(), atomIter1->position());
                  rsq = f.square();
                  type1 = atomIter1->typeId();
                  // Set vector force = (r0-r1)*(forceOverR)
                  if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
                     f *= interactionPtr_->forceOverR(rsq, type0, type1);
                     atomIter0->force() += f;
                     atomIter1->force() -= f;
                  }
               }
            }
         }

         // Iterate over ghosts
         storage().begin(ghostIter);
         if (reverseUpdateFlag()) {

            for ( ; ghostIter.notEnd(); ++ghostIter) {
               id1 = ghostIter->id();
               if (id0 < id1) {
                  if (!atomIter0->mask().isMasked(id1)) {
                     // Set f = r0 - r1, separation between atoms
                     f.subtract(atomIter0->position(), ghostIter->position());
                     rsq = f.square();
                     type1 = ghostIter->typeId();
                     // force = (r0-r1)*(forceOverR)
                     if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
                        f *= interactionPtr_->forceOverR(rsq, type0, type1);
                        atomIter0->force() += f;
                        ghostIter->force() -= f;
                        // Note: If reverseUpdateFlag, increment ghost force 
                     }
                  }
               }
            }

         } else {

            for ( ; ghostIter.notEnd(); ++ghostIter) {
               id1 = ghostIter->id();
               if (!atomIter0->mask().isMasked(id1)) {
                  // Set f = r0 - r1, separation between atoms
                  f.subtract(atomIter0->position(), ghostIter->position());
                  rsq = f.square();
                  type1 = ghostIter->typeId();
                  // force = (r0-r1)*(forceOverR)
                  if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
                     f *= interactionPtr_->forceOverR(rsq, type0, type1);
                     atomIter0->force() += f;
                     // Note: If !reverseUpdateFlag, do not increment ghost force 
                  }
               }
            }

         }

      }
   }

   /*
   * Compute total pair stress (Call on all processors).
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void PairPotentialImpl<Interaction>::computeStress(MPI::Intracomm& communicator)
   #else
   void PairPotentialImpl<Interaction>::computeStress()
   #endif
   {
      // Do nothing if stress is already set.
      if (isStressSet()) return;
 
      Tensor localStress;
      Vector dr;
      Vector f;
      double rsq, forceOverR;
      PairIterator iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int    type0, type1;

      localStress.zero();
      if (reverseUpdateFlag()) {

         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            dr.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = dr.square();
            type0 = atom0Ptr->typeId();
            type1 = atom1Ptr->typeId();
            if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
               f = dr;
               f *= interactionPtr_->forceOverR(rsq, type0, type1);
               incrementPairStress(f, dr, localStress);
            }
         }

      } else {

         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            dr.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = dr.square();
            type0 = atom0Ptr->typeId();
            type1 = atom1Ptr->typeId();
            if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
               f = dr;
               forceOverR = interactionPtr_->forceOverR(rsq, type0, type1);
               assert(!atom0Ptr->isGhost());
               if (atom1Ptr->isGhost()) {
                  forceOverR *= 0.5;
               }
               f *= forceOverR;
               incrementPairStress(f, dr, localStress);
            }
         }

      }

      // Normalize by volume 
      localStress /= boundary().volume();

      // Add localStress from all nodes, set stress to sum on master.
      reduceStress(localStress, communicator);
   }

   /*
   * Compute total pair stress (Call on all processors).
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void PairPotentialImpl<Interaction>::computeForcesAndStress(MPI::Intracomm& communicator)
   #else
   void PairPotentialImpl<Interaction>::computeForcesAndStress()
   #endif
   {
      // If stress is already set, just calculate forces
      if (isStressSet()) {
         computeForces();
         return;
      }
 
      Tensor localStress;
      Vector dr;
      Vector f;
      double rsq;
      PairIterator iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int    type0, type1;

      localStress.zero();
      if (reverseUpdateFlag()) {

         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            dr.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = dr.square();
            type0 = atom0Ptr->typeId();
            type1 = atom1Ptr->typeId();
            if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
               f = dr;
               f *= interactionPtr_->forceOverR(rsq, type0, type1);
               assert(!atom0Ptr->isGhost());
               atom0Ptr->force() += f;
               atom1Ptr->force() -= f;
               incrementPairStress(f, dr, localStress);
            }
         }

      } else {

         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            dr.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = dr.square();
            type0 = atom0Ptr->typeId();
            type1 = atom1Ptr->typeId();
            if (rsq < interactionPtr_->cutoffSq(type0, type1)) {
               f  = dr;
               f *= interactionPtr_->forceOverR(rsq, type0, type1);
               assert(!atom0Ptr->isGhost());
               atom0Ptr->force() += f;
               if (!atom1Ptr->isGhost()) {
                  atom1Ptr->force() -= f;
               } else { // if atom 1 is a ghost
                  f *= 0.5;
               }
               incrementPairStress(f, dr, localStress);
            }
         }

      }

      // Normalize by volume 
      localStress /= boundary().volume();

      // Add localStress from all nodes, set stress to sum on master.
      reduceStress(localStress, communicator);
   }

   /*
   * Compute total pair energies (Call on all processors).
   */
   template <class Interaction>
   #ifdef UTIL_MPI
   void PairPotentialImpl<Interaction>::computePairEnergies(MPI::Intracomm& communicator)
   #else
   void PairPotentialImpl<Interaction>::computePairEnergies()
   #endif
   {
      Vector f;
      double rsq;
      PairIterator iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int    type0, type1;

      DMatrix<double> localPairEnergies;
      localPairEnergies.allocate(nAtomType_, nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         for (int j = 0; j < nAtomType_; ++j) {
            localPairEnergies(i,j) = 0.0;
         }
      }

      //if (methodId() == 0) {
      if (reverseUpdateFlag()) {
         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            assert(!atom0Ptr->isGhost());
            type0 = atom0Ptr->typeId();
            type1 = atom1Ptr->typeId();
            f.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = f.square();
            localPairEnergies(type0, type1) += interactionPtr_->energy(rsq, type0, type1);
         }
      } else {
         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            assert(!atom0Ptr->isGhost());
            type0 = atom0Ptr->typeId();
            type1 = atom1Ptr->typeId();
            f.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = f.square();
            if (!atom1Ptr->isGhost()) {
               localPairEnergies(type0, type1) += interactionPtr_->energy(rsq, type0, type1); 
            } else {
               localPairEnergies(type0, type1) += 0.5*interactionPtr_->energy(rsq, type0, type1); 
            } 
         }
      }

      DMatrix<double> totalPairEnergies;
      totalPairEnergies.allocate(nAtomType_, nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         for (int j = 0; j < nAtomType_; ++j) {
            totalPairEnergies(i,j) = 0.0;
         }
      }

      #ifdef UTIL_MPI
      communicator.Reduce(&localPairEnergies(0,0), &totalPairEnergies(0,0), nAtomType_*nAtomType_,
                           MPI::DOUBLE, MPI::SUM, 0);
      if (communicator.Get_rank() == 0) {
         setPairEnergies(totalPairEnergies);
      } else {
         unsetPairEnergies();
      }
      #else
      setPairEnergies(localPairEnergies);
      #endif
   }

}
#endif
