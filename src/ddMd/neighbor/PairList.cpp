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
    :
      atom1Ptrs_(),
      atom2Ptrs_(),
      first_(),
      cutoff_(0.0),
      atomCapacity_(0),
      pairCapacity_(0),
      nAtom1_(0),
      nAtom2_(0),
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
      int i;

      atomCapacity_ = atomCapacity;
      pairCapacity_ = pairCapacity;
      cutoff_       = cutoff;

      atom1Ptrs_.allocate(atomCapacity_);
      atom2Ptrs_.allocate(pairCapacity_);
      first_.allocate(atomCapacity_ + 1);
   
      // Initialize array elements to null values
      for (i=0; i < atomCapacity_; ++i) {
         atom1Ptrs_[i] = 0;
         first_[i] = PairList::NullIndex;
      }
      first_[atomCapacity_] = PairList::NullIndex;
      for (i=0; i < pairCapacity_; ++i) {
         atom2Ptrs_[i] = 0;
      }

      isAllocated_ = true;
   
   }

   /*
   * Clear the PairList.
   */
   void PairList::clear()
   { 
      nAtom1_   = 0; 
      nAtom2_   = 0; 
      first_[0] = PairList::NullIndex;
   }
 
   /*
   * Build the PairList, i.e., populate it with atom pairs.
   */
   void PairList::build(CellList& cellList, bool reverseUpdateFlag)
   {
      // Precondition
      assert(isAllocated());
 
      Cell::NeighborArray neighbors;
      double dRSq, cutoffSq;
      const Cell*  cellPtr;
      Atom*  atom1Ptr;
      Atom*  atom2Ptr;
      Vector dr;
      int    na;      // number of atoms in this cell
      int    nn;      // number of neighbors for a cell
      int    atom2Id; // atom Id for second atom
      int    i, j;
      bool   foundNeighbor;
  
      // Set maximum squared-separation for pairs in Pairlist
      cutoffSq = cutoff_*cutoff_;
   
      // Initialize counters for primary atoms and neighbors
      nAtom1_   = 0; // Number of primary atoms with neighbors
      nAtom2_   = 0; // Number of secondary atoms (2nd in pair), or pairs
      first_[0] = 0;
   
      // Find all neighbors (cell list)
      cellPtr = cellList.begin();
      while (cellPtr) {
         cellPtr->getNeighbors(neighbors, reverseUpdateFlag);
         na = cellPtr->nAtom();
         nn = neighbors.size();
         for (i = 0; i < na; ++i) {
            atom1Ptr = neighbors[i];
            foundNeighbor = false;

            for (j = 0; j < na; ++j) {
               atom2Ptr = neighbors[j];

               // Neighbors in the same cell as atom1
               if (atom2Ptr > atom1Ptr) {
                  atom2Id  = atom2Ptr->id();
                  if (!atom1Ptr->mask().isMasked(atom2Id)) {
                     dr.subtract(atom2Ptr->position(), atom1Ptr->position()); 
                     dRSq = dr.square();
                     if (dRSq < cutoffSq) {
      
                        if (nAtom2_ >= pairCapacity_) {
                           UTIL_THROW("Overflow: # pairs > pairCapacity_");
                        }
      
                        // If first neighbor of atom1, add atom1 to atom1Ptrs_
                        if (!foundNeighbor) {
   
                           if (nAtom1_ >= atomCapacity_) {
                              UTIL_THROW("Overflow: # pairs > pairCapacity_");
                           }
                           assert(nAtom1_ >= 0);
   
                           atom1Ptrs_[nAtom1_] = atom1Ptr;
                           foundNeighbor = true;
                        }
      
                        // Append 2nd atom to atom2Ptrs_[]
                        atom2Ptrs_[nAtom2_] = atom2Ptr;
                        ++nAtom2_;
      
                     }
                  }
               }
            }

            // Atoms in neighboring cells
            for (j = na; j < nn; ++j) {
               atom2Ptr = neighbors[j];
               atom2Id  = atom2Ptr->id();
               if (!atom1Ptr->mask().isMasked(atom2Id)) {
                  dr.subtract(atom2Ptr->position(), atom1Ptr->position()); 
                  dRSq = dr.square();
                  if (dRSq < cutoffSq) {
   
                     // Check for overflow of atom2Ptrs_
                     if (nAtom2_ >= pairCapacity_) {
                        UTIL_THROW("Overflow: # pairs > pairCapacity_");
                     }
   
                     // If first_ neighbor, record iAtomId in atom1Ptrs_
                     if (!foundNeighbor) {
   
                        if (nAtom1_ >= atomCapacity_) {
                           UTIL_THROW("Overflow: # pairs > pairCapacity_");
                        }
                        assert(nAtom1_ >= 0);
   
                        atom1Ptrs_[nAtom1_] = atom1Ptr;
                        foundNeighbor = true;
                     }
   
                     // Append Id of 2nd atom in pair to atom2Ptrs_[]
                     atom2Ptrs_[nAtom2_] = atom2Ptr;
                     ++nAtom2_;
                 }
              }
   
            }
  
            // When finished with atom1, set next element of first_ array.
            if (foundNeighbor) {
               ++nAtom1_;
               first_[nAtom1_] = nAtom2_;
            }
         }
         cellPtr = cellPtr->nextCellPtr();
      }

      // Increment buildCounter_= number of times the list has been built.
      ++buildCounter_;
 
      // Increment maxima
      if (nAtom1_ > maxNAtomLocal_) maxNAtomLocal_ = nAtom1_;
      if (nAtom2_ > maxNPairLocal_) maxNPairLocal_ = nAtom2_;
     
   }

   /*
   * Initialize a pair iterator.
   */
   void PairList::begin(PairIterator& iterator) const
   {
      iterator.atom1Ptrs_ = &atom1Ptrs_[0];
      iterator.atom2Ptrs_ = &atom2Ptrs_[0];
      iterator.first_     = &first_[0];
      iterator.nAtom1_    = nAtom1_;
      iterator.nAtom2_    = nAtom2_;
      iterator.atom1Id_   = 0;
      iterator.atom2Id_   = 0;
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

   #if 0
   /*
   * Return true if valid, or throw Exception.
   */
   bool PairList::isValid() const
   {
      if (isAllocated()) {

         // Check for null pointers to allocated arrays

         if (atomCapacity_ <= 0) UTIL_THROW("atomCapacity_ <= 0"); 
         if (pairCapacity_ <= 0) UTIL_THROW("pairCapacity_ <= 0"); 

         //if (atom1Ptrs_ == 0)    UTIL_THROW("Null atom1Ptrs_"); 
         //if (atom2Ptrs_ == 0)    UTIL_THROW("Null atom2Ptrs_"); 
         //if (first_ == 0)        UTIL_THROW("Null first_"); 

         if (nAtom1_ < 0) {
            UTIL_THROW("nAtom1_ < 0");
         }
      }

      return true;
   }
   #endif

} 
#endif
