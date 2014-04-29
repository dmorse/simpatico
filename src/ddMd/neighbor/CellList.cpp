#ifndef DDMD_CELL_LIST_CPP
#define DDMD_CELL_LIST_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
      #ifdef UTIL_DEBUG
      maxNAtomCell_(0),
      #endif
      isBuilt_(false)
   {
      for (int i = 0; i < Dimension; ++i) {
         cellLengths_[i] = 0.0;
      }
   }

   /*
   * Destructor. 
   */
   CellList::~CellList()
   {}

   /*
   * Allocate memory for this CellList (generalized coordinates).
   */
   void CellList::allocate(int atomCapacity, const Vector& lower, 
                           const Vector& upper, const Vector& cutoffs, int nCellCut)
   {

      // Allocate arrays of tag and handle objects
      tags_.allocate(atomCapacity);
      atoms_.allocate(atomCapacity);

      // Set grid dimensions and allocate an array of Cell objects
      setGridDimensions(lower, upper, cutoffs, nCellCut);
   }

   /*
   * Calculate number of cells in each direction of grid, resize cells_ array if needed.
   */
   void CellList::setGridDimensions(const Vector& lower, const Vector& upper, 
                                    const Vector& cutoffs, int nCellCut)
   {
      if (nCellCut < 1) {
         UTIL_THROW("Error: nCellCut < 1");
      }
      if (nCellCut > Cell::MaxNCellCut) {
         UTIL_THROW("Error: nCellCut > Cell::MaxNCellCut");
      }
      upper_ = upper;
      lower_ = lower;

      bool isNewGrid;
      if (grid_.size() < 27) {
         isNewGrid = true;
      } else {
         isNewGrid = false;
      }

      Vector lengths;
      IntVector gridDimensions;
      for (int i = 0; i < Dimension; ++i) {
 
         lengths[i] = upper_[i] - lower_[i];
         if (lengths[i] < 0) {
            UTIL_THROW("Processor length[i] < 0.0");
         }
         if (lengths[i] < cutoffs[i]) {
            UTIL_THROW("Processor length[i] < cutoff[i]");
         }
         gridDimensions[i] = int(lengths[i]*nCellCut/cutoffs[i]);
         cellLengths_[i] = lengths[i]/double(gridDimensions[i]);
         lowerOuter_[i] = lower_[i] - nCellCut*cellLengths_[i];
         upperOuter_[i] = upper_[i] + nCellCut*cellLengths_[i];

         // Add two extra layers of cells for ghosts.
         gridDimensions[i] += 2*nCellCut;

         if (gridDimensions[i] != grid_.dimension(i)) {
            isNewGrid = true;   
         }
      }

      // Set new grid dimensions if necessary
      if (isNewGrid) {
         grid_.setDimensions(gridDimensions);
      }

      // Resize and initialize cells_ array, if necessary
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
         // Indicate that cell list must be rebuilt
         isNewGrid = true; 
      }
      assert(newSize >= 27);
      assert(newSize == cells_.size());

      // Build linked cell list, if necessary
      if (isNewGrid) {
         int ic;
         // Initially mark all cells as ghost cells
         for (ic = 0; ic < newSize; ++ic) {
            cells_[ic].setIsGhostCell(true);
         }
         // Loop over local cells, linking and marking each as a local cell.
         IntVector p;
         Cell* prevPtr = 0;
         Cell* cellPtr = 0;
         for (p[0] = nCellCut; p[0] < grid_.dimension(0) - nCellCut; ++p[0]) {
            for (p[1] = nCellCut; p[1] < grid_.dimension(1) - nCellCut; ++p[1]) {
               for (p[2] = nCellCut; p[2] < grid_.dimension(2) - nCellCut; ++p[2]) {
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
      } 
      
   }

   /*
   * Construct grid of cells, build linked list and identify neighbors.
   */
   void CellList::makeGrid(const Vector& lower, const Vector& upper, 
                           const Vector& cutoffs, int nCellCut)
   {

      // Calculate required grid dimensions, reinitialize cells_ array if needed.
      setGridDimensions(lower, upper, cutoffs, nCellCut);

      // Construct e array, to help identify cells within the cutoff.
      // Definition: For i in the range -nCellCut <= i <= nCellCut,
      // let e[i+nCellCut][j] = ( m[i]*celllengths_[j]/cutoffs[j] )**2, 
      // where m[i] = abs(i) - 1 for abs(i) > 0, and m[0] = 0.
      FArray<Vector, 17> e;
      {
         double q, r;
         int i, j;
         for (j = 0; j < Dimension; ++j) {
            q = cellLengths_[j]/cutoffs[j];
            for (i = 1; i <= nCellCut; ++i) {
               r = q*double(i-1);
               r = r*r;
               e[nCellCut + i][j] = r;
               e[nCellCut - i][j] = r;
            }
            e[nCellCut][j] = 0.0;
         }
      }

      // Construct Cell::OffsetArray offsets_ of integer offset strips 
      // Each element strip contains the cell index for the first cell
      // strip.first and the cell index strip.second for the last cell
      // in a contiguous strip of cells for which at least some of the
      // cell lies within a cutoff distance of the primary cell.
      offsets_.clear();
      std::pair<int, int> strip;

      // Add strip (0,0) (self) as the first element of offsets_ array.
      // This guarantees that first nAtom elements in neighborArray are 
      // in the primary cell, allowing for simple self-interaction check.
      strip.first  = 0;
      strip.second = 0;
      offsets_.append(strip); 

      // Loop over all cells within box -nCellCut <= i, j, k <= nCellCut
      double e0, e1, e2;              // Partial sums of distance^2/cutoff^2
      int offset0, offset1, offset;   // Partial sums for cell id offset
      int i, j, k;                    // relative cell coordinates
      const int span0 = grid_.dimension(2)*grid_.dimension(1);
      const int span1 = grid_.dimension(2);
      bool isActive = false; // True iff this cell is within a valid strip
      for (i = -nCellCut; i <= nCellCut; ++i) {
         e0 = e[i+nCellCut][0];
         offset0 = i*span0;
         for (j = -nCellCut; j <= nCellCut; ++j) {
            e1 = e0 + e[j + nCellCut][1];
            offset1 = offset0 + j*span1;
            for (k = -nCellCut; k <= nCellCut; ++k) {
               offset = offset1 + k;
               e2 = e1 + e[k + nCellCut][2];
               if (e2 <= 1.0) {
                  if (offset != 0) { // Exclude offset = 0 (already added)
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
      // Initialize all cells, by associating each with a
      // block of the atoms_ array.

      CellAtom* cellAtomPtr = &atoms_[0];
      for (int i = 0; i < grid_.size(); ++i) {
         cellAtomPtr = cells_[i].initialize(cellAtomPtr);
      }

      // Add all atoms to cells.
      for (int i = 0; i < nAtom_; ++i) {
         cells_[tags_[i].cellRank].append(tags_[i].ptr);
      }

      #ifdef UTIL_DEBUG
      // Calculate maxNAtomCell_
      int nAtomCell;
      maxNAtomCell_ = 0;
      for (int i = 0; i < grid_.size(); ++i) {
         nAtomCell = cells_[i].nAtom();
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
   * Get the maximum number of atoms for which space is allocated.
   */
   int CellList::atomCapacity() const
   {  return tags_.capacity(); }

   /*
   * Get the number of cells for which space has been allocated.
   */
   int CellList::cellCapacity() const
   {  return cells_.capacity(); }

   /*
   * Is this CellList built (i.e., full of atoms)?
   */
   bool CellList::isBuilt() const
   {  return isBuilt_; }

   /*
   * Check validity of CellList, throw an Exception if an error is found.
   */
   bool CellList::isValid() const
   {

      if (isAllocated()) {
         if (tags_.capacity() <= 0) {
            UTIL_THROW("CellList is allocated but tags_.capacity() <= 0");
         }
         if (atoms_.capacity() <= 0) {
            UTIL_THROW("CellList is allocated but atoms_.capacity() <= 0");
         }
      }

      if (isBuilt_) {

         if (!isAllocated()) {
            UTIL_THROW("CellList is built but not allocated");
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
#endif
