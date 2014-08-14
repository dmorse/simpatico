#ifndef SPAN_CELL_LIST_CPP
#define SPAN_CELL_LIST_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CellList.h"
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
//#include <util/containers/FArray.h>

namespace SpAn
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
                             const Vector& upper, const Vector& cutoffs)
   {
      // Allocate arrays of tag and handle objects
      tags_.allocate(atomCapacity);
      atoms_.allocate(atomCapacity);

      // Set grid dimensions and allocate an array of Cell objects
      setGridDimensions(lower, upper, cutoffs);
   }

   /*
   * Calculate number of cells in each direction, resize cells_ array if needed.
   */
   void CellList::setGridDimensions(const Vector& lower, const Vector& upper,
                                      const Vector& cutoffs)
   {
      upper_ = upper;
      lower_ = lower;

      bool isNewGrid;
      if (grid_.size() < 27) {
         isNewGrid = true;
      } else {
         isNewGrid = false;
      }

      // Calculate gridDimensions = number of cells in each direction
      // If any differ from current grid dimension, isNewGrid = false.
      Vector lengths;
      IntVector gridDimensions;
      for (int i = 0; i < Dimension; ++i) {
         lengths[i] = upper_[i] - lower_[i];
         if (lengths[i] < 0) {
            UTIL_THROW("Error: length[i] < 0.0");
         }
         if (lengths[i] < cutoffs[i]) {
            UTIL_THROW("Error: length[i] < cutoff[i]");
         }
         gridDimensions[i] = int(lengths[i]/cutoffs[i]);
         cellLengths_[i] = lengths[i]/double(gridDimensions[i]);

         if (gridDimensions[i] != grid_.dimension(i)) {
            isNewGrid = true;
         }
      }

      // Set new grid dimensions if necessary
      if (isNewGrid) {
         grid_.setDimensions(gridDimensions);

         gridOffsets_[Dimension - 1] = 1;
         for (int i = Dimension - 1; i > 0; --i) {
            gridOffsets_[i-1] = gridOffsets_[i]*gridDimensions[i];
         }
      }

      // Resize and initialize cells_ array, if necessary
      int oldSize = cells_.size();
      int newSize = grid_.size();
      if (newSize != oldSize) {
         cells_.resize(newSize);
         if (newSize > oldSize) {
            for (int i = 0; i < newSize; ++i) {
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
         // Loop over local cells, linking cells.
         IntVector p;
         Cell* prevPtr = 0;
         Cell* cellPtr = 0;
         int ic;
         for (p[0] = 0; p[0] < grid_.dimension(0); ++p[0]) {
            for (p[1] = 0; p[1] < grid_.dimension(1); ++p[1]) {
               for (p[2] = 0; p[2] < grid_.dimension(2); ++p[2]) {
                  ic = grid_.rank(p);
                  cellPtr = &cells_[ic];
                  if (prevPtr) {
                     prevPtr->setNextCell(*cellPtr);
                  }
                  else {
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
                             const Vector& cutoffs)
   {
      // Calculate required grid dimensions, reinitialize cells_ array if needed.
      setGridDimensions(lower, upper, cutoffs);

      // Calculate offsets to move to neighboring cells
      for (int cellId = 0; cellId < grid_.size(); cellId++) {
         Cell::OffsetArray *offsets = new Cell::OffsetArray;

         for (int i = 0; i < Dimension; i++) {
            std::pair<int, int> axisOffset;
            axisOffset.first  = calculateAxisOffset(cellId, i,  1);
            axisOffset.second = calculateAxisOffset(cellId, i, -1);

            offsets->append(axisOffset);
         }

         cells_[cellId].setOffsetArray(*offsets);
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
   * Get maximum number of atoms in any cell. 
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
