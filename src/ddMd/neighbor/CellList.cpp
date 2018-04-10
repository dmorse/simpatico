/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CellList.h"
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/containers/FArray.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   CellList::CellList()
    : begin_(0),
      nAtom_(0),
      nReject_(0),
      nCellCut_(-1),
      #ifdef UTIL_DEBUG
      maxNAtomCell_(0),
      #endif
      isBuilt_(false)
   {
      IntVector gridDimensions;
      for (int i = 0; i < Dimension; ++i) {
         cellLengths_[i] = 0.0;
         gridDimensions[i] = 1;
      }
      grid_.setDimensions(gridDimensions);
      UTIL_CHECK(grid_.size() == 1);
   }

   /*
   * Destructor.
   */
   CellList::~CellList()
   {}

   /*
   * Allocate or reallocate memory for atoms in this CellList.
   */
   void CellList::setAtomCapacity(int atomCapacity)
   {
      if (tags_.capacity() == 0) {
         UTIL_CHECK(atoms_.capacity() == 0);
         tags_.allocate(atomCapacity);
         atoms_.allocate(atomCapacity);
      } else {
         UTIL_CHECK(tags_.capacity() == atoms_.capacity());
         if (atomCapacity > tags_.capacity()) {
            tags_.deallocate();
            atoms_.deallocate();
            tags_.allocate(atomCapacity);
            atoms_.allocate(atomCapacity);
         }
      }
   }

   /*
   * Calculate number of cells in each direction of grid, resize
   * cells_ array and reconstruct linked list if the grid changes,
   * and recompute the offsets_ array.
   */
   void CellList::makeGrid(const Vector& lower, const Vector& upper,
                           const Vector& cutoffs, int nCellCut)
   {
      if (nCellCut < 1) {
         UTIL_THROW("Error: nCellCut < 1");
      }
      if (nCellCut > Cell::MaxNCellCut) {
         UTIL_THROW("Error: nCellCut > Cell::MaxNCellCut");
      }

      // Set default value of isNewGrid
      bool isNewGrid = false;
      if (grid_.size() < 27) {
         isNewGrid = true;
      } else 
      if (nCellCut_ == -1) {
         isNewGrid = true;
      } else
      if (nCellCut_ != nCellCut) {
         isNewGrid = true;
      }

      // Set member variables to values given as parameters
      upper_ = upper;
      lower_ = lower;
      nCellCut_ = nCellCut;

      // Compute required grid dimensions.
      Vector lengths;            // dimensions of region owned by processor
      IntVector gridDimensions;  // grid dimensions, including ghost regions 
      for (int i = 0; i < Dimension; ++i) {

         lengths[i] = upper_[i] - lower_[i];
         if (lengths[i] < 0) {
            UTIL_THROW("Processor length[i] < 0.0");
         }
         if (lengths[i] < cutoffs[i]) {
            // Note: This would cause parallelization strategy to fail.
            UTIL_THROW("Processor length[i] < cutoff[i]");
         }
         gridDimensions[i] = int(lengths[i]*nCellCut_/cutoffs[i]);
         cellLengths_[i] = lengths[i]/double(gridDimensions[i]);

         // Add two extra layers of cells for ghosts.
         gridDimensions[i] += 2*nCellCut_;
         lowerOuter_[i] = lower_[i] - nCellCut_*cellLengths_[i];
         upperOuter_[i] = upper_[i] + nCellCut_*cellLengths_[i];

         if (gridDimensions[i] != grid_.dimension(i)) {
            isNewGrid = true;
         }
      }

      /*
      * Here, isNewGrid is true iff any grid dimension has changed or
      * nCellCut has changed. 
      */

      if (isNewGrid) {

         // Set new dimensions in grid_ object
         grid_.setDimensions(gridDimensions);

         // Resize and initialize cells_ array
         int oldSize = cells_.size();
         int newSize = grid_.size();
         if (newSize != oldSize) {
            cells_.resize(newSize);
            if (newSize > oldSize) {
               for (int i = 0; i < newSize; ++i) {
                  cells_[i].setOffsetArray(offsets_);
                  cells_[i].setId(i);
               }
            }
         }

         // Initially mark all cells as ghost cells
         int ic; // cell index
         for (ic = 0; ic < newSize; ++ic) {
            cells_[ic].setIsGhostCell(true);
         }

         // Create linked list of local cells, and mark each as local.
         IntVector p;                        // grid coordinate indices
         IntVector d = grid_.dimensions();   // grid dimensions
         Cell* prevPtr = 0;
         Cell* cellPtr = 0;
         for (p[0] = nCellCut_; p[0] < d[0] - nCellCut_; ++p[0]) {
            for (p[1] = nCellCut_; p[1] < d[1] - nCellCut_; ++p[1]) {
               for (p[2] = nCellCut_; p[2] < d[2] - nCellCut_; ++p[2]) {
                  ic = grid_.rank(p);
                  cellPtr = &cells_[ic];
                  cellPtr->setIsGhostCell(false);
                  if (prevPtr) {
                     prevPtr->setNextCell(*cellPtr);
                  } else {
                     begin_ = cellPtr;
                  }
                  prevPtr = cellPtr;
               }
            }
         }
         cellPtr->setLastCell();

      } // end if(isNewGrid)

      /*
      * The remainder of this function computes the offsets_ array.
      * The offsets_ array is recomputed whenever makeGrid is called
      * because the offsets_ array can change due to changes in the 
      * cutoffs vector even when the grid dimensions and nCellCut are 
      * unchanged. When coordinates are expressed in scaled form, for 
      * which each coordinate in the primary unit cell is in the range 
      * [0,1], changes in physical unit cell size or shape during an
      * NPT simulations cause changes in the elements of the elements 
      * of the Vector cutoffs. In simulations with a rigid unit cell,
      * the makeGrid function should only be called at the beginning
      * of the simulation, rather than every step.
      */

      /*
      * Construct e array, to help identify cells within the cutoff.
      * Definition: For i in the range -nCellCut_ <= i <= nCellCut_, 
      * let e[i+nCellCut_][j] = ( m[i]*celllengths_[j]/cutoffs[j] )**2,
      * where m[i] = abs(i) - 1 for abs(i) > 0, and m[0] = 0.
      */
      FArray<Vector, 17> e;
      {
         double q, r;
         int i, j;
         for (j = 0; j < Dimension; ++j) {
            q = cellLengths_[j]/cutoffs[j];
            for (i = 1; i <= nCellCut_; ++i) {
               r = q*double(i-1);
               r = r*r;
               e[nCellCut_ + i][j] = r;
               e[nCellCut_ - i][j] = r;
            }
            e[nCellCut_][j] = 0.0;
         }
      }
  
      /* 
      * Construct Cell::OffsetArray offsets_ of integer offset strips.
      * Each element of the offsets_ array is a pair of integers 
      * representing relative indices for a contiguous strip of cells 
      * for which at least some of the cell lies within a cutoff 
      * distance of the primary cell. The first integer in each pair 
      * is the relative index of the first cell in such a strip, and 
      * the second is the relative index of the last cell in the strip.
      */
      offsets_.clear();
      std::pair<int, int> strip;
  
      /* 
      * Add pair (0,0) (self) as the first element of offsets_ array.
      * This guarantees that first nAtom elements in neighborArray are
      * in the primary cell, allowing for simple self-interaction check.
      */
      strip.first  = 0;
      strip.second = 0;
      offsets_.append(strip);
   
      // Loop over all cells within box -nCellCut_ <= i, j, k <= nCellCut_
      double e0, e1, e2;             // Partial sums of distance^2/cutoff^2
      int offset0, offset1, offset;  // Partial sums for cell id offset
      int i, j, k;                   // relative cell coordinates
      const int span0 = grid_.dimension(2)*grid_.dimension(1);
      const int span1 = grid_.dimension(2);
      bool isActive = false; // True iff cell is within a valid strip
      for (i = -nCellCut_; i <= nCellCut_; ++i) {
         e0 = e[i+nCellCut_][0];
         offset0 = i*span0;
         for (j = -nCellCut_; j <= nCellCut_; ++j) {
            e1 = e0 + e[j + nCellCut_][1];
            offset1 = offset0 + j*span1;
            for (k = -nCellCut_; k <= nCellCut_; ++k) {
               offset = offset1 + k;
               e2 = e1 + e[k + nCellCut_][2];
               if (e2 <= 1.0) {
                  if (offset != 0) { 
                     // Note: Exclude offset = 0, because already added
                     if (isActive) {
                        if (offset == strip.second + 1) {
                           strip.second = offset;
                        } else {
                           offsets_.append(strip);
                           strip.first  = offset;
                           strip.second = offset;
                        }
                     } else {
                        strip.first  = offset;
                        strip.second = offset;
                        isActive = true;
                     }
                  } else {
                     if (isActive) {
                        offsets_.append(strip);
                        isActive = false;
                     }
                  }
               } else {
                  if (isActive) {
                     offsets_.append(strip);
                     isActive = false;
                  }
               }
            }
         }
      }

      // Append last strip to offsets_, if still active at end of loop.
      if (isActive) {
         offsets_.append(strip);
      }

   }

   /*
   * Resets all cells to empty state.
   */
   void CellList::clear()
   {
      // Clear all Cell objects
      if (grid_.size() > 0) {
         for (int i = 0; i < grid_.size(); ++i) {
            cells_[i].clear();
         }
      }

      nAtom_ = 0;
      nReject_ = 0;
      #ifdef UTIL_DEBUG
      maxNAtomCell_ = 0;
      #endif

      isBuilt_ = false;
   }

   /*
   * Initialize all cells and fill them with atoms.
   */
   void CellList::build()
   {
      // Associate each cell with a block of the atoms_ array.
      CellAtom* cellAtomPtr = &atoms_[0];
      for (int i = 0; i < grid_.size(); ++i) {
         cellAtomPtr = cells_[i].initialize(cellAtomPtr);
      }

      // Add all atoms to cells.
      for (int i = 0; i < nAtom_; ++i) {
         cells_[tags_[i].cellRank].append(tags_[i].ptr);
      }

      // Note: Cell::append() calls CellAtom::setPtr() to set the pointer
      // to the appropriate atom, but does not call CellAtom::update(). 

      #ifdef UTIL_DEBUG
      // Calculate maxNAtomCell_
      int nAtomCell;
      maxNAtomCell_ = 0;
      for (int i = 0; i < grid_.size(); ++i) {
         nAtomCell = cells_[i].nAtom();
         UTIL_CHECK(nAtomCell == cells_[i].atomCapacity());
         if (nAtomCell > maxNAtomCell_) {
            maxNAtomCell_ = nAtomCell;
         }
      }
      #endif

      isBuilt_ = true;
   }

   /*
   * Update position information in all CellAtom objects.
   */
   void CellList::update()
   {
      for (int i = 0; i < nAtom_; ++i) {
         atoms_[i].update();
      }
   }

   /*
   * Get total number of atoms in this CellList.
   */
   int CellList::nAtom() const
   {  return nAtom_; }

   /*
   * Get number of atoms that were rejected (not added to cells).
   */
   int CellList::nReject() const
   {  return nReject_; }

   #ifdef UTIL_DEBUG
   /*
   * Get number of atoms that were rejected (not added to cells).
   */
   int CellList::maxNAtomCell() const
   {  return maxNAtomCell_; }
   #endif

   /*
   * Check validity of CellList, throw an Exception if an error is found.
   */
   bool CellList::isValid() const
   {

      // Array dimensions sanity checks
      UTIL_CHECK(cells_.capacity() >= 0);
      UTIL_CHECK(atoms_.capacity() >= 0);
      UTIL_CHECK(cells_.capacity() >= 0);
      UTIL_CHECK(atoms_.capacity() == tags_.capacity());

      // If grid has been built
      if (grid_.size() > 1) {
         UTIL_CHECK(grid_.size() >= 27);
         UTIL_CHECK(cells_.capacity() >= grid_.size());
         UTIL_CHECK(nCellCut_ >= 1);
      }

      if (isBuilt_) {

         if (cells_.capacity() <= 0) {
            UTIL_THROW("CellList is built but cells_ is not allocated");
         }
         if (atoms_.capacity() <= 0) {
            UTIL_THROW("CellList is built but atoms_ is not allocated");
         }

         // Check validity of all cells individually.
         // const CellAtom* atomPtr;
         const Cell* cellPtr;
         int   nAtomCell;
         int   nAtomSum = 0;
         for (int icell = 0; icell < grid_.size(); ++icell) {
            cellPtr = &cells_[icell];
            nAtomCell = cellPtr->nAtom();
            if (nAtomCell != cellPtr->atomCapacity()) {
               UTIL_THROW("Cell nAtom != atomCapacity");
            }
            #if 0
            if (nAtomCell > 0) {
               for (int i = 0; i < nAtomCell; ++i) {
                  atomPtr = cellPtr->atomPtr(i);
                  if (icell != cellIndexFromPosition(atomPtr->position())) {
                      UTIL_THROW("Inconsistent position");
                  }
               }
            }
            #endif
            nAtomSum += nAtomCell;
         }

         // Check that total number of atoms in all cells equals nAtom.
         // Note: nAtom is incremented by the placeAtom() method.
         if (nAtom_ >= 0) {
            if (nAtomSum != nAtom_) {
               UTIL_THROW("Number of atoms in all cells != nAtom");
            }
         }

      } else { // if not isBuilt_

         if (nAtom_ != 0) {
            UTIL_THROW("CellList is not built, but nAtom_ != 0");
         }

      }

      return true;
   }

}
