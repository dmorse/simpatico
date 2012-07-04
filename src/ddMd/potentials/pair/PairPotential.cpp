#ifndef DDMD_PAIR_POTENTIAL_CPP
#define DDMD_PAIR_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "PairPotential.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/neighbor/PairIterator.h>
#include <ddMd/neighbor/PairList.h>
#include <ddMd/neighbor/CellList.h>
#include <ddMd/communicate/Domain.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace DdMd
{
   using namespace Util;

   /*
   * Default constructor (for unit testing).
   */
   PairPotential::PairPotential()
    : skin_(0.0),
      cutoff_(0.0),
      pairCapacity_(0),
      domainPtr_(0),
      boundaryPtr_(0),
      storagePtr_(0),
      timer_(PairPotential::NTime),
      methodId_(0),
      nPair_(0),
      forceCommFlag_(false)
   {} 

   /*
   * Constructor.
   */
   PairPotential::PairPotential(Simulation& simulation)
    : skin_(0.0),
      cutoff_(0.0),
      pairCapacity_(0),
      domainPtr_(&simulation.domain()),
      boundaryPtr_(&simulation.boundary()),
      storagePtr_(&simulation.atomStorage()),
      timer_(PairPotential::NTime),
      methodId_(0),
      nPair_(0),
      forceCommFlag_(false)
   {}

   /*
   * Associate with related objects. (for unit testing).
   */
   void PairPotential::associate(Domain& domain, Boundary& boundary, 
                                 AtomStorage& storage)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;
      storagePtr_ = &storage;
   } 

   /*
   * Destructor.
   */
   PairPotential::~PairPotential()
   {}

   /*
   * Set flag to specify if reverse force communication is enabled.
   */
   void PairPotential::setForceCommFlag(bool forceCommFlag)
   {  forceCommFlag_ = forceCommFlag; }

   void PairPotential::readPairListParam(std::istream& in)
   {
      read<double>(in, "skin", skin_);
      read<int>(in, "pairCapacity", pairCapacity_);
      read<Boundary>(in, "maxBoundary", maxBoundary_);

      // Set upper and lower bound of the processor domain.
      boundary().setLengths(maxBoundary_.lengths());
      Vector lower;
      Vector upper;
      for (int i = 0; i < Dimension; ++i) {
         lower[i] = domain().domainBound(i, 0);
         upper[i] = domain().domainBound(i, 1);
      }

      // Allocate CellList and PairList
      int atomCapacity = storage().atomCapacity() + storage().atomCapacity();
      cutoff_ = maxPairCutoff() + skin_;
      cellList_.allocate(atomCapacity, lower, upper, cutoff_);
      pairList_.allocate(atomCapacity, pairCapacity_, cutoff_);
   }

   /*
   * Allocate memory for the cell list.
   */
   void PairPotential::initialize(const Vector& lower, const Vector& upper, 
                                double skin, int pairCapacity)
   {
      skin_ = skin;
      pairCapacity_ = pairCapacity;

      cutoff_ = maxPairCutoff() + skin;
      int atomCapacity = storage().atomCapacity();

      cellList_.allocate(atomCapacity, lower, upper, cutoff_);
      pairList_.allocate(atomCapacity, pairCapacity_, cutoff_);
   }

   /*
   * Build the cell list (i.e., fill with atoms).
   */
   void PairPotential::findNeighbors(const Vector& lower, const Vector& upper)
   {
      stamp(PairPotential::START);

      cellList_.makeGrid(lower, upper, cutoff_);
      cellList_.clear();
     
      // Add all atoms to the cell list. 
      AtomIterator atomIter;
      storage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         cellList_.placeAtom(*atomIter);
      }

      // Add all ghosts to the cell list. 
      GhostIterator ghostIter;
      storage().begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         cellList_.placeAtom(*ghostIter);
      }

      cellList_.build();
      assert(cellList_.isValid());
      assert(cellList_.nAtom() + cellList_.nReject() == storage().nAtom() + storage().nGhost());
      stamp(PairPotential::BUILD_CELL_LIST);

      pairList_.build(cellList_, forceCommFlag());
      stamp(PairPotential::BUILD_PAIR_LIST);
   }

   /*
   * Build the cell list (i.e., fill with atoms).
   */
   void PairPotential::findNeighbors()
   {
      // Make the cell list grid.
      Vector lower;
      Vector upper;
      for (int i = 0; i < Dimension; ++i) {
         lower[i] = domain().domainBound(i, 0);
         upper[i] = domain().domainBound(i, 1);
      }
      findNeighbors(lower, upper);
   }

   /*
   * Compute total pair nPair on all processors.
   */
   #ifdef UTIL_MPI
   void PairPotential::computeNPair(MPI::Intracomm& communicator)
   #else
   void PairPotential::computeNPair()
   #endif
   { 
      double cut = maxPairCutoff();
      double cutSq = cut*cut;
      int localNPair = 0; 
      if (methodId() == 0) {
         localNPair = nPairList(cutSq); 
      } else 
      if (methodId() == 1) {
         localNPair = nPairCell(cutSq); 
      } else {
         localNPair = nPairNSq(cutSq); 
      }
      nPair_ = localNPair;
      #ifdef UTIL_MPI
      communicator.Reduce(&localNPair, &nPair_, 1, MPI::INT, MPI::SUM, 0);
      #else
      nPair_ = localNPair;
      #endif
   }

   /**
   * Return twice the number of pairs within the specified cutoff.
   * 
   * This method should only be called on the rank 0 processor. The
   * return value is computed by a previous call to computeNPair.
   */
   int PairPotential::nPair() const
   {  return nPair_; }

   /*
   * Count number pairs using pair list (private).
   */
   int PairPotential::nPairList(double cutoffSq)
   {
      Vector f;
      double rsq;
      PairIterator iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int count = 0; // atom counter
      if (forceCommFlag()) {
         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            assert(!atom0Ptr->isGhost());
            f.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = f.square();
            if (rsq < cutoffSq) {
               count += 2;
            }
         }
      } else {
         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            assert(!atom0Ptr->isGhost());
            f.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = f.square();
            if (rsq < cutoffSq) {
               if (!atom1Ptr->isGhost()) {
                  count += 2;
               } else {
                  count += 1;
               }
            }
         }
      }
      return count;
   }

   /*
   * Count number pairs using cell list (private).
   */
   int PairPotential::nPairCell(double cutoffSq)
   {
      // Find all neighbors (cell list)
      Cell::NeighborArray neighbors;
      Vector f;
      double rsq;
      Atom*  atomPtr0;
      Atom*  atomPtr1;
      const Cell*  cellPtr;
      int na, nn, i, j;
      int count = 0;

      // Iterate over linked list of local cells.
      cellPtr = cellList_.begin();
      while (cellPtr) {
         cellPtr->getNeighbors(neighbors, forceCommFlag());
         na = cellPtr->nAtom();
         nn = neighbors.size();
         for (i = 0; i < na; ++i) {
            atomPtr0 = neighbors[i];

            // Loop over atoms in this cell
            for (j = 0; j < na; ++j) {
               atomPtr1 = neighbors[j];
               if (atomPtr1 > atomPtr0) {
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  if (rsq < cutoffSq) {
                     count += 2;
                  }
               }
            }

            // Loop over atoms in neighboring cells.
            if (forceCommFlag()) {
               for (j = na; j < nn; ++j) {
                  atomPtr1 = neighbors[j];
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  if (rsq < cutoffSq) {
                     count += 2;
                  }
               }
            } else {
               for (j = na; j < nn; ++j) {
                  atomPtr1 = neighbors[j];
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  if (rsq < cutoffSq) {
                     if (!atomPtr1->isGhost()) {
                        count += 2;
                     } else {
                        count += 1;
                     }
                  }
               }
            }

         }
         cellPtr = cellPtr->nextCellPtr();
      } // while (cellPtr) 
      return count;
   }

   /*
   * Count pairs using n-squared loop. (private)
   */
   int PairPotential::nPairNSq(double cutoffSq)
   {
      Vector f;
      double rsq;
      AtomIterator  atomIter0, atomIter1;
      GhostIterator ghostIter;
      int id0, id1;
      int count = 0;

      // Iterate over local atom 0
      storage().begin(atomIter0);
      for ( ; atomIter0.notEnd(); ++atomIter0) {
         id0 = atomIter0->id();

         // Iterate over local atom 1
         storage().begin(atomIter1);
         for ( ; atomIter1.notEnd(); ++atomIter1) {
            id1 = atomIter1->id();
            if (id0 < id1) {
               if (!atomIter0->mask().isMasked(id1)) {
                  f.subtract(atomIter0->position(), atomIter1->position());
                  rsq = f.square();
                  if (rsq < cutoffSq) {
                     count += 2;
                  }
               }
            }
         }

         // Iterate over ghost atoms
         storage().begin(ghostIter);
         if (forceCommFlag()) {
            for ( ; ghostIter.notEnd(); ++ghostIter) {
               id1 = ghostIter->id();
               if (id0 < id1) {
                  if (!atomIter0->mask().isMasked(id1)) {
                     f.subtract(atomIter0->position(), ghostIter->position());
                     rsq = f.square();
                     if (rsq < cutoffSq) {
                        count += 2;
                     }
                  }
               }
            }
         } else {
            for ( ; ghostIter.notEnd(); ++ghostIter) {
               id1 = ghostIter->id();
               if (!atomIter0->mask().isMasked(id1)) {
                  f.subtract(atomIter0->position(), ghostIter->position());
                  rsq = f.square();
                  if (rsq < cutoffSq) {
                     count += 1;
                  }
               }
            }
         }

      }
      return count;
   }

}
#endif
