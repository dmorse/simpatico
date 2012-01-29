#ifndef INTERACTION_CPP
#define INTERACTION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Interaction.h"
#include <ddMd/system/System.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/neighbor/PairIterator.h>
#include <ddMd/communicate/Domain.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   Interaction::Interaction()
    : skin_(0.0),
      cutoff_(0.0),
      storagePtr_(0),
      potentialPtr_(0),
      domainPtr_(0),
      boundaryPtr_(0),
      pairCapacity_(0),
      methodId_(0)
   {}

   /*
   * Constructor.
   */
   Interaction::Interaction(System& system)
    : skin_(0.0),
      cutoff_(0.0),
      storagePtr_(&system.atomStorage()),
      potentialPtr_(&system.pairPotential()),
      domainPtr_(&system.domain()),
      boundaryPtr_(&system.boundary()),
      pairCapacity_(0),
      methodId_(0)
   {}

   /*
   * Destructor.
   */
   Interaction::~Interaction()
   {}

   /*
   * Retain pointers to associated objects.
   */
   void Interaction::associate(AtomStorage& storage, const PairPotential& potential)
   { 
      storagePtr_   = &storage;
      potentialPtr_ = &potential;
   }

   void Interaction::readParam(std::istream& in)
   {
      readBegin(in,"Interaction");
      read<double>(in, "skin", skin_);
      read<int>(in, "pairCapacity", pairCapacity_);
      read<Boundary>(in, "maxBoundary", maxBoundary_);

      cutoff_ = potentialPtr_->maxPairCutoff() + skin_;
      int atomCapacity = storagePtr_->atomCapacity();
                       + storagePtr_->ghostCapacity();

      // Set upper and lower bound of this domain.
      boundaryPtr_->setLengths(maxBoundary_.lengths());
      Vector lower;
      Vector upper;
      for (int i = 0; i < Dimension; ++i) {
         lower[i] = domainPtr_->domainBound(i, 0);
         upper[i] = domainPtr_->domainBound(i, 1);
      }
      cellList_.allocate(atomCapacity, lower, upper, cutoff_);
      pairList_.allocate(atomCapacity, pairCapacity_, cutoff_);

      readEnd(in);
   }

   /*
   * Allocate memory for the cell list.
   */
   void Interaction::setParam(const Vector& lower, const Vector& upper, 
                              double skin, int pairCapacity)
   {
      skin_ = skin;
      pairCapacity_ = pairCapacity;

      cutoff_ = potentialPtr_->maxPairCutoff() + skin;
      int atomCapacity = storagePtr_->atomCapacity();
                       + storagePtr_->ghostCapacity();

      cellList_.allocate(atomCapacity, lower, upper, cutoff_);
      pairList_.allocate(atomCapacity, pairCapacity_, cutoff_);
   }

   void Interaction::setMethodId(int methodId)
   {  methodId_ = methodId; }

   /*
   * Build the cell list (i.e., fill with atoms).
   */
   void Interaction::findNeighbors(const Vector& lower, const Vector& upper)
   {
      cellList_.makeGrid(lower, upper, cutoff_);
      cellList_.clear();
     
      // Add all atoms to the cell list. 
      AtomIterator atomIter;
      storagePtr_->begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         cellList_.placeAtom(*atomIter);
      }

      // Add all ghosts to the cell list. 
      GhostIterator ghostIter;
      storagePtr_->begin(ghostIter);
      for ( ; !ghostIter.atEnd(); ++ghostIter) {
         cellList_.placeAtom(*ghostIter);
      }

      cellList_.build();
      assert(cellList_.isValid());

      pairList_.build(cellList_);
   }

   /*
   * Build the cell list (i.e., fill with atoms).
   */
   void Interaction::findNeighbors()
   {
      // Make the cell list grid.
      Vector lower;
      Vector upper;
      for (int i = 0; i < Dimension; ++i) {
         lower[i] = domainPtr_->domainBound(i, 0);
         upper[i] = domainPtr_->domainBound(i, 1);
      }
      findNeighbors(lower, upper);
   }

   /*
   * Set forces on all local atoms to zero.
   */
   void Interaction::zeroForces()
   {
      AtomIterator atomIter;
      storagePtr_->begin(atomIter); 
      for( ; !atomIter.atEnd(); ++atomIter){
         atomIter->force().zero();
      }
   }

   /*
   * Set forces on all local atoms to zero.
   */
   void Interaction::calculateForces()
   {
      zeroForces();
      addPairForces();
   }

   /*
   * Increment atomic forces, without calculating energy.
   */
   void Interaction::addPairForces()
   {  
       if (methodId_ == 0) {
          addPairForcesList(true, false); 
       } else
       if (methodId_ == 1) {
          addPairForcesCell(true, false); 
       } else {
          addPairForcesNSq(true, false); 
       }
   }

   /*
   * Increment atomic forces and compute pair energy.
   */
   void Interaction::addPairForces(double& energy)
   {  
       if (methodId_ == 0) {
          energy = addPairForcesList(true, true); 
       } else
       if (methodId_ == 1) {
          energy = addPairForcesCell(true, true); 
       } else {
          energy = addPairForcesNSq(true, true); 
       }
   }

   /*
   * Calculate pair energy.
   */
   double Interaction::pairPotential()
   {  
       if (methodId_ == 0) {
          return addPairForcesList(false, true); 
       } else 
       if (methodId_ == 1) {
          return addPairForcesCell(false, true); 
       } else {
          return addPairForcesNSq(false, true); 
       }
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   double Interaction::addPairForcesList(bool needForce, bool needEnergy)
   {
      Vector f;
      double rsq, energy;
      PairIterator iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int    type0, type1;

      for (pairList_.begin(iter); !iter.atEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         type0 = atom0Ptr->typeId();
         type1 = atom1Ptr->typeId();
         f.subtract(atom0Ptr->position(), atom1Ptr->position());
         rsq = f.square();
         if (!atom1Ptr->isGhost()) {
            if (needEnergy) {
               energy += potentialPtr_->energy(rsq, type0, type1);
            }
            if (needForce) {
               f *= potentialPtr_->forceOverR(rsq, type0, type1);
               atom0Ptr->force() += f;
               atom1Ptr->force() -= f;
            }
         } else {
            if (needEnergy) {
               energy += 0.5*potentialPtr_->energy(rsq, type0, type1);
            }
            if (needForce) {
               f *= potentialPtr_->forceOverR(rsq, type0, type1);
               atom0Ptr->force() += f;
            }
         }
      }
      return energy;
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   double Interaction::addPairForcesCell(bool needForce, bool needEnergy)
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
                  if (needEnergy) {
                     energy += potentialPtr_->energy(rsq, type0, type1);
                  }
                  if (needForce) {
                     f *= potentialPtr_->forceOverR(rsq, type0, type1);
                     atomPtr0->force() += f;
                     atomPtr1->force() -= f;
                  }
               }
            }

            // Loop over atoms in neighboring cells.
            for (j = na; j < nn; ++j) {
               atomPtr1 = neighbors[j];
               type1 = atomPtr1->typeId();
               f.subtract(atomPtr0->position(), atomPtr1->position());
               rsq = f.square();
               if (!atomPtr1->isGhost()) {
                  if (needEnergy) {
                     energy += potentialPtr_->energy(rsq, type0, type1);
                  }
                  if (needForce) {
                     f *= potentialPtr_->forceOverR(rsq, type0, type1);
                     atomPtr0->force() += f;
                     atomPtr1->force() -= f;
                  }
               } else {
                  if (needEnergy) {
                     energy += 0.5*potentialPtr_->energy(rsq, type0, type1);
                  }
                  if (needForce) {
                     f *= potentialPtr_->forceOverR(rsq, type0, type1);
                     atomPtr0->force() += f;
                  }
               }
            }

         }

         cellPtr = cellPtr->nextCellPtr();

      } // while (cellPtr) 

      return energy;
   }

   /*
   * Increment atomic forces and/or pair energy (private).
   */
   double Interaction::addPairForcesNSq(bool needForce, bool needEnergy)
   {
      // Preconditions
      //if (!storagePtr_->isInitialized()) {
      //   UTIL_THROW("AtomStorage must be initialized");
      //}

      Vector f;
      double rsq;
      double energy = 0.0;
      AtomIterator  atomIter0, atomIter1;
      GhostIterator ghostIter;
      int           type0, type1;

      // Iterate over atom 0
      storagePtr_->begin(atomIter0);
      for ( ; !atomIter0.atEnd(); ++atomIter0) {
         type0 = atomIter0->typeId();

         // Iterate over atom 1
         storagePtr_->begin(atomIter1);
         for ( ; !atomIter1.atEnd(); ++atomIter1) {
            type1 = atomIter1->typeId();

            if (atomIter0->id() < atomIter1->id()) {

               // Set f = r0 - r1, separation between atoms
               f.subtract(atomIter0->position(), atomIter1->position());
               rsq = f.square();
 
               if (needEnergy) {
                  energy += potentialPtr_->energy(rsq, type0, type1);
               }

               if (needForce) {
 
                  // Set vector force = (r0-r1)*(forceOverR)
                  f *= potentialPtr_->forceOverR(rsq, type0, type1);
           
                  // Add equal and opposite forces.
                  atomIter0->force() += f;
                  atomIter1->force() -= f;

               }

            }
            
         }

         // Iterate over ghosts
         storagePtr_->begin(ghostIter);
         for ( ; !ghostIter.atEnd(); ++ghostIter) {
            type1 = ghostIter->typeId();

            // Set f = r0 - r1, separation between atoms
            f.subtract(atomIter0->position(), ghostIter->position());
            rsq = f.square();

            // Set force = (r0-r1)*(forceOverR)
            f *= potentialPtr_->forceOverR(rsq, type0, type1);
     
            // Add half energy of local-ghost interaction. 
            if (needEnergy) {
               energy += 0.5*potentialPtr_->energy(rsq, type0, type1);
            }
 
            // Add force to local atom
            atomIter0->force() += f;
            
         }

      }

      return energy;
   }

}
#endif
