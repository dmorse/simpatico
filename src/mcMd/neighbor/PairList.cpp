#ifndef MCMD_PAIR_LIST_CPP
#define MCMD_PAIR_LIST_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "PairList.h"
#include "PairIterator.h"
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace McMd
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
      oldPositions_(),
      skin_(0.0),
      cutoff_(0.0),
      atomCapacity_(0),
      pairCapacity_(0),
      nAtom1_(0),
      nAtom2_(0),
      nAtom_(0),
      tList1_(0),
      maxNAtom_(0),
      maxNAtom2_(0),
      buildCounter_(0)
   {  setClassName("PairList"); }
   
   /*
   * Destructor.
   */
   PairList::~PairList() 
   {}

   /*
   * Read atomCapacity and pairCapacity from file.
   */
   void PairList::readParameters(std::istream& in) 
   {
      read<int>(in, "atomCapacity", atomCapacity_);
      read<int>(in, "pairCapacity", pairCapacity_);
      read<double>(in, "skin", skin_);
   }

   /*
   * Load parameters an archive, but do not allocate.
   */
   void PairList::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar, "atomCapacity", atomCapacity_);
      loadParameter<int>(ar, "pairCapacity", pairCapacity_);
      loadParameter<double>(ar, "skin", skin_);
      ar >> cutoff_;
      ar >> cellList_;

      // Allocate neighbor list data structures
      atom1Ptrs_.allocate(atomCapacity_);
      atom2Ptrs_.allocate(pairCapacity_);
      first_.allocate(atomCapacity_ + 1);
      oldPositions_.allocate(atomCapacity_);
   
      // Initialize array elements to null values
      for (int i=0; i < atomCapacity_; ++i) {
         atom1Ptrs_[i] = 0;
         first_[i] = PairList::NullIndex;
         oldPositions_[i].zero();
      }
      first_[atomCapacity_] = PairList::NullIndex;
      for (int i=0; i < pairCapacity_; ++i) {
         atom2Ptrs_[i] = 0;
      }
      clear();
   }

   /*
   * Save parameters an archive.
   */
   void PairList::save(Serializable::OArchive &ar)
   {
      ar << atomCapacity_;
      ar << pairCapacity_;
      ar << skin_;
      ar << cutoff_;
      ar << cellList_;
   }

   /*
   * Allocate CellList and PairList arrays, initialize to empty state.
   */
   void PairList::allocate(int atomIdEnd, const Boundary& boundary, 
                                          double potentialCutoff)
   {
      int i;

      // Preconditions
      if (skin_ < 0.0000001) {
         UTIL_THROW("skin must be set before PairList::allocate");
      }
      if (atomCapacity_ <= 0 ) {
         UTIL_THROW("atomCapacity_ must be set before PairList::allocate");
      }
      if (pairCapacity_ <= 0 ) {
         UTIL_THROW("pairCapacity_ must be set before PairList::allocate");
      }

      // PairList cutoff = cutoff for potential + a "skin"
      cutoff_ = potentialCutoff + skin_;

      // Initialize the private CellList
      cellList_.allocate(atomIdEnd, boundary, cutoff_);
   
      // Allocate neighbor list data structures
      atom1Ptrs_.allocate(atomCapacity_);
      atom2Ptrs_.allocate(pairCapacity_);
      first_.allocate(atomCapacity_ + 1);
      oldPositions_.allocate(atomCapacity_);
   
      // Initialize array elements to null values
      for (i=0; i < atomCapacity_; ++i) {
         atom1Ptrs_[i] = 0;
         first_[i] = PairList::NullIndex;
         oldPositions_[i].zero();
      }
      first_[atomCapacity_] = PairList::NullIndex;
      for (i=0; i < pairCapacity_; ++i) {
         atom2Ptrs_[i] = 0;
      }
   }

   /*
   * Make the grid of cells for the internal CellList.
   */
   void PairList::makeGrid(const Boundary& boundary)
   {  cellList_.makeGrid(boundary, cutoff_); }
 
   /*
   * Clear the CellList and PairList.
   */
   void PairList::clear()
   { 
      cellList_.clear(); 
      nAtom1_   = 0; 
      nAtom2_   = 0; 
      nAtom_    = 0;
      tList1_   = atomCapacity_ - 1;
      first_[0] = 0;
   }
 
   /*
   * Build the PairList, i.e., populate it with atom pairs.
   */
   void PairList::build(const Boundary& boundary)
   {
      CellList::NeighborArray cellNeighbor;
      Vector  iPos, jPos;
      Atom   *iAtomPtr, *jAtomPtr;
      double  dRSq, cutoffSq;
      int     nCellNeighbor, nCellAtom, totCells;
      int     ic, ip, iAtomId, jp, jAtomId;
      bool    foundNeighbor;
  
      // Precondition
      assert(isAllocated());
 
      // Set maximum squared-separation for pairs in Pairlist
      cutoffSq = cutoff_*cutoff_;
   
      // Initialize counters for primary atoms and neighbors
      nAtom1_   = 0; // Number of primary atoms with neighbors
      nAtom2_   = 0; // Number of secondary atoms (2nd in pair), or pairs
      nAtom_    = 0; // Number of atoms, with or without neighbors
      tList1_   = atomCapacity_ - 1;
      first_[0] = 0;
   
      // Loop over cells containing primary atom. ic = cell index
      totCells = cellList_.totCells();
      for (ic = 0; ic < totCells; ++ic) {
   
         // Get Array cellNeighbor of Ids of neighbor atoms for cell ic.
         // Elements 0,..., nCellAtom - 1 contain Ids for atoms in cell ic.
         // Elements nCellAtom,..., nCellNeighbor-1 are from neighboring cells.
         cellList_.getCellNeighbors(ic, cellNeighbor, nCellAtom);
         nCellNeighbor = cellNeighbor.size();
      
         // Loop over atoms in cell ic
         for (ip = 0; ip < nCellAtom; ++ip) {
            iAtomPtr = cellNeighbor[ip]; 
            iPos     = iAtomPtr->position();
            iAtomId  = iAtomPtr->id();
            ++nAtom_;
            if (nAtom_ > atomCapacity_) {
               UTIL_THROW("Overflow: nAtom_ > atomCapacity_ in PairList");
            }
   
            // Indicate that no neighbors have been found yet
            foundNeighbor = false;
   
            // Loop over atoms in all neighboring cells, including cell ic.
            for (jp = 0; jp < nCellNeighbor; ++jp) {
               jAtomPtr = cellNeighbor[jp]; 
               jPos     = jAtomPtr->position();
               jAtomId  = jAtomPtr->id();
        
               // Avoid double counting: only count pairs with jAtomId > iAtomId
               if ( jAtomId > iAtomId ) {
     
                  // Exclude bonded pairs
                  if (!iAtomPtr->mask().isMasked(*jAtomPtr))  {
   
                     // Calculate distance between atoms i and j
                     dRSq = boundary.distanceSq(iPos, jPos);
   
                     if (dRSq < cutoffSq) {
      
                        // Check for overflow of atom2Ptrs_
                        if (nAtom2_ >= pairCapacity_) {
                           UTIL_THROW("Overflow: # pairs > pairCapacity_");
                        }
      
                        // If first_ neighbor, record iAtomId in atom1Ptrs_
                        if (!foundNeighbor) {

                           assert(nAtom1_ < atomCapacity_);
                           assert(nAtom1_ >= 0);
                           assert(tList1_ >= nAtom1_);

                           atom1Ptrs_[nAtom1_] = iAtomPtr;
                           oldPositions_[nAtom1_] = iPos;
                           foundNeighbor = true;
                        }
      
                        // Append Id of 2nd atom in pair to atom2Ptrs_[]
                        assert(tList1_ <  atomCapacity_);
                        atom2Ptrs_[nAtom2_] = jAtomPtr;
                        ++nAtom2_;
      
                     }
                  }
   
               } // end if jAtomId > iAtomId
         
            } // end for jp (j atom)
   
            if (foundNeighbor) {
               // Update counter first_[nAtom1_]
               ++nAtom1_;
               first_[nAtom1_] = nAtom2_;
            } else {
               // Add atom with no neighbors to top of atom1Ptrs_ 
               assert(tList1_ >= nAtom1_);
               assert(tList1_ <  atomCapacity_);
               assert(tList1_ >= 0);
               atom1Ptrs_[tList1_] = iAtomPtr;
               oldPositions_[tList1_] = iPos;
               --tList1_;
            }
         
         } // end for ip (i atom)
   
      } // end for ic (i cell)

      // Increment buildCounter_= number of times the list has been built.
      ++buildCounter_;
 
      // Increment maximum
      if (nAtom_  > maxNAtom_)  maxNAtom_  = nAtom_;
      if (nAtom2_ > maxNAtom2_) maxNAtom2_ = nAtom2_;
     
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
   * Return false if any atom has moved a distance greater than skin/2,
   * or if the PairList has never been built. Return true otherwise.
   */
   bool PairList::isCurrent(const Boundary& boundary) const
   {
      // Precondition
      assert(isAllocated());
 
      double dRSq, dRSqMax;
      int    ip;
  
      // If the list has never been built, it is not current. 
      if (buildCounter_ == 0) return false;
   
      dRSqMax = 0.25*skin_*skin_;

      // Loop over atoms with neighbors (increment from bottom)
      for (ip = 0; ip < nAtom1_; ++ip) {
         dRSq = boundary.distanceSq(atom1Ptrs_[ip]->position(), 
                                    oldPositions_[ip]);
         if (dRSq > dRSqMax) {
            return false;
         }
      }
   
      // Loop over atoms with no neighbors (decrement from top)
      for (ip = atomCapacity_ - 1; ip > tList1_; --ip) {
         dRSq = boundary.distanceSq(atom1Ptrs_[ip]->position(), 
                                    oldPositions_[ip]);
         if (dRSq > dRSqMax) {
            return false;
         }
      }

      return true;
   }

   /*
   * Clear all statistics.
   */
   void PairList::clearStatistics()
   {
      maxNAtom_     = 0;
      maxNAtom2_    = 0;
      buildCounter_ = 0;
   }


   /*
   * Return true if valid, or throw Exception.
   */
   bool PairList::isValid() const
   {
      if (isAllocated()) {

         // Check for null pointers to allocated arrays

         if (atomCapacity_ <= 0) UTIL_THROW("atomCapacity_ <= 0"); 
         if (pairCapacity_ <= 0) UTIL_THROW("atomCapacity_ <= 0"); 

         //if (atom1Ptrs_ == 0)    UTIL_THROW("Null atom1Ptrs_"); 
         //if (atom2Ptrs_ == 0)    UTIL_THROW("Null atom2Ptrs_"); 
         //if (first_ == 0)        UTIL_THROW("Null first_"); 
         //if (oldPositions_ == 0) UTIL_THROW("Null oldPositions_"); 

         if (nAtom_ > 0) {
            if (nAtom_ != nAtom1_ + atomCapacity_ - 1 -  tList1_) {
               UTIL_THROW("Inconsistent innAtom_, nAtom1_, tList1");
            }
            cellList_.isValid(nAtom_);
            if (nAtom_ != cellList_.nAtom()) {
               UTIL_THROW("Inconsistent values of nAtom_");
            }
         } else 
         if (nAtom_ == 0) {
            if (nAtom1_ != 0) UTIL_THROW("nAtom_ ==0 and nAtom1 != 0");
            if (nAtom2_ != 0) UTIL_THROW("nAtom_ ==0 and nAtom2 != 0");
            cellList_.isValid();
         } else {
            UTIL_THROW("nAtom_ < 0");
         }
      }

      return true; 
   }

} 
#endif
