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
#include <util/mpi/MpiLoader.h>
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
      methodId_(0),
      nPair_(0),
      pairEnergies_()
   {  setClassName("PairPotential"); } 

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
      methodId_(0),
      nPair_(0),
      pairEnergies_()
   {  setClassName("PairPotential"); } 

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
   * Allocate memory for the cell list.
   */
   void PairPotential::initialize(const Boundary& maxBoundary, double skin, 
                                  int pairCapacity)
   {
      skin_ = skin;
      pairCapacity_ = pairCapacity;
      maxBoundary_ = maxBoundary;
      cutoff_ = maxPairCutoff() + skin;

      allocate();
   }

   /*
   * Read parameters for PairList and allocate memory.  
   */
   void PairPotential::readParameters(std::istream& in)
   {
      read<double>(in, "skin", skin_);
      read<int>(in, "pairCapacity", pairCapacity_);
      read<Boundary>(in, "maxBoundary", maxBoundary_);
      cutoff_ = maxPairCutoff() + skin_;
      allocate();
   }

   /*
   * Read parameters for PairList and allocate memory.  
   */
   void PairPotential::loadParameters(Serializable::IArchive& ar)
   {
  
      loadParameter<double>(ar, "skin", skin_);
      loadParameter<int>(ar, "pairCapacity", pairCapacity_);
      loadParameter<Boundary>(ar, "maxBoundary", maxBoundary_);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(cutoff_);
      loader.load(methodId_);
      allocate();
   }

   /*
   * Save parameters to output/saving Archive.
   */
   void PairPotential::save(Serializable::OArchive& ar)
   {
      ar << skin_;
      ar << pairCapacity_;
      ar << maxBoundary_;
      ar << cutoff_;
      ar << methodId_;
   }

   /*
   * Allocate memory for the cell list and pair list.
   *
   * On entry, cutoff_, pairCapacity_, and maxBoundary_ must be defined.
   */
   void PairPotential::allocate()
   {
      // Allocate PairList
      int localCapacity = storage().atomCapacity();
      pairList_.allocate(localCapacity, pairCapacity_, cutoff_);

      // Calculate cell list cutoff lengths for all directions
      Vector cutoffs;
      Vector lower;
      Vector upper;
      for (int i = 0; i < Dimension; ++i) {
         lower[i] = domain().domainBound(i, 0);
         upper[i] = domain().domainBound(i, 1);
         cutoffs[i] = cutoff_/maxBoundary_.length(i);
      }

      // Allocate CellList
      int totalCapacity = localCapacity + storage().ghostCapacity();
      cellList_.allocate(totalCapacity, lower, upper, cutoffs);
   }

   /*
   * Build the cell list.
   */
   void PairPotential::buildCellList()
   {
      if (storage().isCartesian()) {
         UTIL_THROW("Coordinates are Cartesian entering buildCellList");
      }

      // Set cutoff and domain bounds.
      Vector cutoffs;
      Vector lower;
      Vector upper;
      for (int i = 0; i < Dimension; ++i) {
         cutoffs[i] = cutoff_/boundaryPtr_->length(i);
         lower[i] = domain().domainBound(i, 0);
         upper[i] = domain().domainBound(i, 1);
      }
      cellList_.makeGrid(lower, upper, cutoffs);
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

      // Finalize cell list
      cellList_.build();
     
      // Postconditions
      assert(cellList_.isValid());
      assert(cellList_.nAtom() + cellList_.nReject() 
             == storage().nAtom() + storage().nGhost());
      if (storage().isCartesian()) {
        UTIL_THROW("Coordinates are Cartesian exiting buildCellList");
      }
   }

   /*
   * Build the pair list.
   */
   void PairPotential::buildPairList()
   {
      if (!storage().isCartesian()) {
         UTIL_THROW("Coordinates not Cartesian entering buildPairList");
      }
      pairList_.build(cellList_, reverseUpdateFlag());
   }

   /*
   * Return value of pair energies.
   */
   DMatrix<double> PairPotential::pairEnergies() const
   {  return pairEnergies_.value(); }

   /*
   * Set a value for pair energies.
   */
   void PairPotential::setPairEnergies(DMatrix<double> pairEnergies)
   {  pairEnergies_.set(pairEnergies); }

   /*
   * Mark pair energies as unknown. 
   */
   void PairPotential::unsetPairEnergies()
   {  pairEnergies_.unset(); }

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

   /*
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
      if (reverseUpdateFlag()) {
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
         cellPtr->getNeighbors(neighbors, reverseUpdateFlag());
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
            if (reverseUpdateFlag()) {
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
         if (reverseUpdateFlag()) {
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
