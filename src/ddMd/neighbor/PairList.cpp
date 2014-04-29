#ifndef DDMD_PAIR_LIST_CPP
#define DDMD_PAIR_LIST_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "PairList.h"
#include "PairIterator.h"
#include <ddMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   PairList::PairList()
    : atom1Ptrs_(),
      atom2Ptrs_(),
      first_(),
      cutoff_(0.0),
      atomCapacity_(0),
      pairCapacity_(0),
      maxNAtomLocal_(0),
      maxNPairLocal_(0),
      buildCounter_(0),
      maxNAtom_(0),
      maxNPair_(0),
      isAllocated_(false)
   {}
   
   /*
   * Destructor.
   */
   PairList::~PairList() 
   {}
   
   /*
   * Allocate CellList and PairList arrays, initialize to empty state.
   */
   void PairList::allocate(int atomCapacity, int pairCapacity , double cutoff) 
   {
      atomCapacity_ = atomCapacity;
      pairCapacity_ = pairCapacity;
      cutoff_       = cutoff;

      atom1Ptrs_.reserve(atomCapacity_);
      atom2Ptrs_.reserve(pairCapacity_);
      first_.reserve(atomCapacity_ + 1);
  
      isAllocated_ = true;
   }

   /*
   * Clear the PairList.
   */
   void PairList::clear()
   { 
      atom1Ptrs_.clear();
      atom2Ptrs_.clear();
      first_.clear();
   }
 
   /*
   * Build the PairList, i.e., populate it with atom pairs.
   */
   void PairList::build(CellList& cellList, bool reverseUpdateFlag)
   {
      // Precondition
      assert(isAllocated());
 
      Cell::NeighborArray neighbors;
      double cutoffSq;
      const Cell* cellPtr;
      CellAtom* atom1Ptr;
      CellAtom* atom2Ptr;
      Mask* maskPtr;
      Vector dr;
      int na;                 // number of atoms in this cell
      int nn;                 // number of neighbors for a cell
      int i, j;
      bool hasNeighbor;
  
      // Set maximum squared-separation for pairs in Pairlist
      cutoffSq = cutoff_*cutoff_;
   
      // Initialize counters for primary atoms and neighbors
      atom1Ptrs_.clear();
      atom2Ptrs_.clear();
      first_.clear();
      first_.append(0);

      // Copy positions and ids into cell list
      cellList.update();
   
      // Find all neighbors (cell list)
      cellPtr = cellList.begin();
      while (cellPtr) {
         na = cellPtr->nAtom(); // # of atoms in cell

         if (na) {
            cellPtr->getNeighbors(neighbors, reverseUpdateFlag);
            nn = neighbors.size();
   
            // Loop over primary atoms (atom1) in primary cell
            for (i = 0; i < na; ++i) {
               atom1Ptr = neighbors[i];
               maskPtr  = atom1Ptr->maskPtr();
   
               // Loop over secondary atoms
               hasNeighbor = false;
               for (j = i + 1; j < nn; ++j) {
                  atom2Ptr = neighbors[j];
                  dr.subtract(atom2Ptr->position(), atom1Ptr->position()); 
                  if (dr.square() < cutoffSq && !maskPtr->isMasked(atom2Ptr->id())) {
                     atom2Ptrs_.append(atom2Ptr->ptr());
                     hasNeighbor = true;
                  }
               }

               // Complete processing of atom1.
               if (hasNeighbor) {
                  atom1Ptrs_.append(atom1Ptr->ptr());
                  first_.append(atom2Ptrs_.size());
               }

            } // for ia 
         } // if (na)

         // Advance to next cell in a linked list
         cellPtr = cellPtr->nextCellPtr();
      }

      // Postconditions
      if (atom1Ptrs_.size()) {
         if (first_.size() != atom1Ptrs_.size() + 1) {
            UTIL_THROW("Array size problem");
         }
         if (first_[0] != 0) {
            UTIL_THROW("Incorrect first element of first_");
         }
         if (first_[atom1Ptrs_.size()] != atom2Ptrs_.size()) {
            UTIL_THROW("Incorrect last element of first_");
         }
      }

      // Increment buildCounter_= number of times the list has been built.
      ++buildCounter_;
 
      // Increment maxima
      if (atom1Ptrs_.size() > maxNAtomLocal_) {
         maxNAtomLocal_ = atom1Ptrs_.size();
      }
      if (atom2Ptrs_.size() > maxNPairLocal_) {
         maxNPairLocal_ = atom2Ptrs_.size();
      }
   }

   /*
   * Initialize a pair iterator.
   */
   void PairList::begin(PairIterator& iterator) const
   {
      if (atom1Ptrs_.size()) {
         iterator.atom1Ptrs_ = &atom1Ptrs_[0];
         iterator.atom2Ptrs_ = &atom2Ptrs_[0];
         iterator.first_     = &first_[0];
         iterator.nAtom1_    = atom1Ptrs_.size();
         iterator.nAtom2_    = atom2Ptrs_.size();
         iterator.atom1Id_   = 0;
         iterator.atom2Id_   = 0;
      }
   }

   /*
   * Compute memory usage statistics (call on all processors).
   */
   #ifdef UTIL_MPI
   void PairList::computeStatistics(MPI::Intracomm& communicator)
   #else
   void PairList::computeStatistics()
   #endif
   { 
      #ifdef UTIL_MPI
      int maxNAtomGlobal;
      int maxNPairGlobal;
      communicator.Allreduce(&maxNAtomLocal_, &maxNAtomGlobal, 1, 
                             MPI::INT, MPI::MAX);
      communicator.Allreduce(&maxNPairLocal_, &maxNPairGlobal, 1, 
                             MPI::INT, MPI::MAX);
      maxNAtom_.set(maxNAtomGlobal);
      maxNPair_.set(maxNPairGlobal);
      maxNAtomLocal_ = maxNAtomGlobal;
      maxNPairLocal_ = maxNPairGlobal;
      #else
      maxNAtom_.set(maxNAtomLocal_);
      maxNPair_.set(maxNPairLocal_);
      #endif
   }

   /*
   * Clear all statistics.
   */
   void PairList::clearStatistics() 
   {
      maxNAtomLocal_ = 0;
      maxNPairLocal_ = 0;
      maxNAtom_.unset();
      maxNPair_.unset();
      buildCounter_ = 0;
   }

   /*
   * Output statistics.
   */
   void PairList::outputStatistics(std::ostream& out)
   {

      out << std::endl;
      out << "PairList" << std::endl;
      out << "NPair: max, capacity     " 
                  << Int(maxNPair_.value(), 10)
                  << Int(pairCapacity_, 10)
                  << std::endl;
      out << "NAtom: max, capacity     " 
                  << Int(maxNAtom_.value(), 10)
                  << Int(atomCapacity_, 10)
                  << std::endl;
   }

} 
#endif
