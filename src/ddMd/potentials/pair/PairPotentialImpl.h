#ifndef DDMD_PAIR_POTENTIAL_IMPL_H
#define DDMD_PAIR_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/potentials/pair/PairPotential.h>
#include <util/global.h>

namespace Util {
   class Vector;
   class Tensor;
}

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
      * Set integer id to specify algorithm.
      *
      * \param methodId algorithm id: 0=pair list, 1=cell list, 2=N^2 loop.
      */
      void setMethodId(int methodId);
  
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
      virtual void readParam(std::istream& in);

      /// \name Interaction interface
      //@{

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNAtomType(int nAtomType);
  
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
      virtual void addForces();

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
      * Get the total energy calculated previously by computeEnergy().
      *
      * Call only on master. 
      */
      virtual double energy();

      #if 0 
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
      #endif

      //@}

   private:

      /**
      * Total pair energy on all processors (valid only on master).
      */
      double energy_;
 
      /**
      * Pointer to pair interaction object.
      */ 
      Interaction* interactionPtr_;

      // Index for method used to calculate forces / energies.
      int methodId_;

      /**
      * Calculate atomic pair energy, using PairList.
      */
      double energyList();

      /**
      * Calculate atomic pair forces, using PairList.
      */
      void addForcesList();

      /**
      * Calculate atomic pair forces and/or pair potential energy.
      */
      double energyCell();

      /**
      * Calculate atomic pair forces and/or pair potential energy.
      */
      void addForcesCell();

      /**
      * Calculate atomic pair energy, using N^2 loop.
      * 
      * Use an O(N^2) double loop over all atoms.
      */
      double energyNSq();

      /**
      * Calculate atomic pair forces, using N^2 loop.
      */
      void addForcesNSq();

      #if 0 
      template <typename T>
      void computeStressImpl(T& stress) const;
      #endif

   };

}

#include <ddMd/simulation/Simulation.h>
//#include <ddMd/simulation/stress.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/neighbor/PairIterator.h>
#include <ddMd/communicate/Domain.h>

#include <util/space/Dimension.h>
//#include <util/space/Vector.h>
#include <util/space/Tensor.h>
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
      interactionPtr_(0),
      methodId_(0)
   {  interactionPtr_ = new Interaction; }
 
   /* 
   * Default constructor.
   */
   template <class Interaction>
   PairPotentialImpl<Interaction>::PairPotentialImpl()
    : PairPotential(),
      interactionPtr_(0),
      methodId_(0)
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

   /*
   * Set parameter to determine which method to use to calculate forces.
   */
   template <class Interaction>
   void PairPotentialImpl<Interaction>::setMethodId(int methodId)
   {  methodId_ = methodId; }

   template <class Interaction>
   void PairPotentialImpl<Interaction>::readParam(std::istream& in)
   {
      readBegin(in,"PairPotential");
      bool nextIndent = false;
      readParamComposite(in, interaction(), nextIndent);
      readPairListParam(in);
      readEnd(in);
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
   PairPotentialImpl<Interaction>::energy(double rsq, 
                                        int iAtomType, int jAtomType) const
   { return interaction().energy(rsq, iAtomType, jAtomType); }

   /*
   * Return force / separation for a single pair.
   */
   template <class Interaction> double 
   PairPotentialImpl<Interaction>::forceOverR(double rsq, 
                                        int iAtomType, int jAtomType) const
   { return interaction().forceOverR(rsq, iAtomType, jAtomType); }

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
   void PairPotentialImpl<Interaction>::addForces()
   {  
       stamp(PairPotential::START);
       if (methodId_ == 0) {
          addForcesList(); 
       } else
       if (methodId_ == 1) {
          addForcesCell(); 
       } else {
          addForcesNSq(); 
       }
       stamp(PairPotential::FORCES);
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
      double localEnergy = 0; 
      if (methodId_ == 0) {
         localEnergy = energyList(); 
      } else 
      if (methodId_ == 1) {
         localEnergy = energyCell(); 
      } else {
         localEnergy = energyNSq(); 
      }

      #ifdef UTIL_MPI
      communicator.Reduce(&localEnergy, &energy_, 1, 
                          MPI::DOUBLE, MPI::SUM, 0);
      #else
      energy_ = localEnergy;
      #endif
   }

   /*
   * Return total pair energy from all processors.
   */
   template <class Interaction>
   double PairPotentialImpl<Interaction>::energy()
   {  return energy_; } 

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
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
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
      return energy;
   }


   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction>
   void PairPotentialImpl<Interaction>::addForcesList()
   {
      Vector f;
      double rsq;
      PairIterator iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int    type0, type1;
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         f.subtract(atom0Ptr->position(), atom1Ptr->position());
         rsq = f.square();
         type0 = atom0Ptr->typeId();
         type1 = atom1Ptr->typeId();
         f *= interactionPtr_->forceOverR(rsq, type0, type1);
         atom0Ptr->force() += f;
         if (!atom1Ptr->isGhost()) {
            atom1Ptr->force() -= f;
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
                  energy += interactionPtr_->energy(rsq, type0, type1);
               }
            }
            // Loop over atoms in neighboring cells.
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
         cellPtr = cellPtr->nextCellPtr();
      } // while (cellPtr) 
      return energy;
   }

   /*
   * Increment atomic forces using Cell List (private).
   */
   template <class Interaction>
   void PairPotentialImpl<Interaction>::addForcesCell()
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
                  f *= interactionPtr_->forceOverR(rsq, type0, type1);
                  atomPtr0->force() += f;
                  atomPtr1->force() -= f;
               }
            }
            // Loop over atoms in neighboring cells.
            for (j = na; j < nn; ++j) {
               atomPtr1 = neighbors[j];
               type1 = atomPtr1->typeId();
               f.subtract(atomPtr0->position(), atomPtr1->position());
               rsq = f.square();
               f *= interactionPtr_->forceOverR(rsq, type0, type1);
               atomPtr0->force() += f;
               if (!atomPtr1->isGhost()) {
                  atomPtr1->force() -= f;
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
      int           type0, type1;

      // Iterate over local atom 0
      storage().begin(atomIter0);
      for ( ; atomIter0.notEnd(); ++atomIter0) {
         type0 = atomIter0->typeId();
         // Iterate over local atom 1
         storage().begin(atomIter1);
         for ( ; atomIter1.notEnd(); ++atomIter1) {
            type1 = atomIter1->typeId();
            if (atomIter0->id() < atomIter1->id()) {
               f.subtract(atomIter0->position(), atomIter1->position());
               rsq = f.square();
               energy += interactionPtr_->energy(rsq, type0, type1);
            }
         }
         // Iterate over ghost atoms
         storage().begin(ghostIter);
         for ( ; ghostIter.notEnd(); ++ghostIter) {
            type1 = ghostIter->typeId();
            f.subtract(atomIter0->position(), ghostIter->position());
            rsq = f.square();
            energy += 0.5*interactionPtr_->energy(rsq, type0, type1);
         }
      }
      return energy;
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   template <class Interaction>
   void PairPotentialImpl<Interaction>::addForcesNSq()
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

         // Iterate over atom 1
         storage().begin(atomIter1);
         for ( ; atomIter1.notEnd(); ++atomIter1) {
            type1 = atomIter1->typeId();
            id1 =   atomIter1->id();
            if (id0 < id1) {
               if (!atomIter0->mask().isMasked(id1)) {
                  // Set f = r0 - r1, separation between atoms
                  f.subtract(atomIter0->position(), atomIter1->position());
                  rsq = f.square();
                  // Set vector force = (r0-r1)*(forceOverR)
                  f *= interactionPtr_->forceOverR(rsq, type0, type1);
                  atomIter0->force() += f;
                  atomIter1->force() -= f;
               }
            }
         }

         // Iterate over ghosts
         storage().begin(ghostIter);
         for ( ; ghostIter.notEnd(); ++ghostIter) {
            id1 = ghostIter->id();
            type1 = ghostIter->typeId();
            if (!atomIter0->mask().isMasked(id1)) {
               // Set f = r0 - r1, separation between atoms
               f.subtract(atomIter0->position(), ghostIter->position());
               rsq = f.square();
               // force = (r0-r1)*(forceOverR)
               f *= interactionPtr_->forceOverR(rsq, type0, type1);
               atomIter0->force() += f;
            }
         }

      }
   }

   #if 0
   /* 
   * Add nonBonded pair forces to atomic forces.
   */
   template <class Interaction>
   template <typename T>
   void PairPotentialImpl<Interaction>::computeStressImpl(T& stress) const
   {

      Vector       dr;
      Vector       force;
      double       rsq;
      PairIterator iter;
      Atom        *atom1Ptr;
      Atom        *atom0Ptr;

      setToZero(stress);

      // Loop over nonbonded neighbor pairs
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position(), dr);
         force  = dr;
         force *= interaction().
                  forceOverR(rsq, atom0Ptr->typeId(), atom1Ptr->typeId());
         incrementPairStress(force, dr, stress);
      }

      // Normalize by volume 
      stress /= boundary().volume();
      normalizeStress(stress);

   }

   template <class Interaction>
   void PairPotentialImpl<Interaction>::computeStress(double& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Interaction>
   void PairPotentialImpl<Interaction>::computeStress(Util::Vector& stress)
        const
   {  computeStressImpl(stress); }

   template <class Interaction>
   void PairPotentialImpl<Interaction>::computeStress(Util::Tensor& stress)
        const
   {  computeStressImpl(stress); }
   #endif

}
#endif
