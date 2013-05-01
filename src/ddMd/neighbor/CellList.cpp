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
    : nAtom_(0),
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
      handles_.allocate(atomCapacity);

      // Allocate array of Cell objects
      setGridDimensions(lower, upper, cutoffs);
      cells_.allocate(grid_.size());

      for (int i = 0; i < grid_.size(); ++i) {
         cells_[i].setOffsetArray(offsets_);
         cells_[i].setId(i);
      }

   }

   /*
   * Allocate memory for this CellList.
   */
   void CellList::allocate(int atomCapacity, const Vector& lower, 
                           const Vector& upper, double cutoff)
   {
      Vector cutoffs;
      for (int i = 0; i < Dimension; ++i) {
         cutoffs[i] = cutoff;
      }
      allocate(atomCapacity, lower, upper, cutoffs);
   }

   /*
   * Calculate number of cells in each direction of grid.
   */
   void CellList::setGridDimensions(const Vector& lower, const Vector& upper, 
                                    const Vector& cutoffs)
   {
      Vector    lengths;
      IntVector gridDimensions;
      upper_ = upper;
      lower_ = lower;

      for (int i = 0; i < Dimension; ++i) {

         lengths[i] = upper_[i] - lower_[i];
         assert(lengths[i] > 0.0);
         gridDimensions[i] = int(lengths[i]/cutoffs[i]);
         if (gridDimensions[i] < 1) {
            gridDimensions[i] = 1;
         }
         cellLengths_[i] = lengths[i]/double(gridDimensions[i]);
         lowerOuter_[i] = lower_[i] - cellLengths_[i];
         upperOuter_[i] = upper_[i] + cellLengths_[i];

         // Add two extra layers of cells for ghosts.
         gridDimensions[i] += 2;
   
      }
      grid_.setDimensions(gridDimensions);

      // Postcondition
      if (grid_.size() < 1) {
         UTIL_THROW("totCells_ must be >= 1");
      }

   }

   /*
   * Construct grid of cells, build linked list and identify neighbors.
   */
   void CellList::makeGrid(const Vector& lower, const Vector& upper, const Vector& cutoffs)
   {

      // Calculate required grid dimensions
      setGridDimensions(lower, upper, cutoffs);

      if (grid_.size() > cells_.capacity()) {
         UTIL_THROW("Insufficient memory was allocated for this grid");
      }

      // Mark all cells to as ghost cells by default.
      int ic;
      for (ic = 0; ic < grid_.size(); ++ic) {
         cells_[ic].setIsGhostCell(true);
      }

      // Build linked list of local cells.
      IntVector p;
      Cell* prevPtr = 0;
      Cell* cellPtr = 0;
      for (p[0] = 1; p[0] < grid_.dimension(0) - 1; ++p[0]) {
         for (p[1] = 1; p[1] < grid_.dimension(1) - 1; ++p[1]) {
            for (p[2] = 1; p[2] < grid_.dimension(2) - 1; ++p[2]) {
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

      // Calculate range of displacements to neighbor cells
      IntVector dmin, dmax;
      for (int i = 0; i < Dimension; ++i) {
         if (grid_.dimension(i) > 2) {
            dmin[i] = -1;
            dmax[i] =  1;
         } else if (grid_.dimension(i) == 2) {
            dmin[i] =  0;
            dmax[i] =  1;
         } else if (grid_.dimension(i) == 1) {
            dmin[i] =  0;
            dmax[i] =  0;
         }
      }

      // Construct array of integer offsets to neighbors
      IntVector span;
      IntVector d;
      int       offset;
      span[2] = 1;
      span[1] = grid_.dimension(2);
      span[0] = grid_.dimension(2)*grid_.dimension(1);
      offsets_.clear();
      offsets_.append(0);
      for (d[0] = dmin[0]; d[0] <= dmax[0]; ++d[0]) {
         for (d[1] = dmin[1]; d[1] <= dmax[1]; ++d[1]) {
            for (d[2] = dmin[2]; d[2] <= dmax[2]; ++d[2]) {
               offset = d[0]*span[0] + d[1]*span[1] + d[2];
               if (offset != 0) {
                  offsets_.append(offset);
               }
            }
         }
      }

   }

   /*
   * Construct grid of cells, build linked list and identify neighbors (Cartesian).
   */
   void CellList::makeGrid(const Vector& lower, const Vector& upper, double cutoff)
   {
      Vector cutoffs;
      for (int i = 0; i < Dimension; ++i) {
         cutoffs[i] = cutoff;
      }
      makeGrid(lower, upper, cutoffs);
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
      // block of the handles_ array.

      Atom** handlePtr = &handles_[0];
      for (int i = 0; i < grid_.size(); ++i) {
         handlePtr = cells_[i].initialize(handlePtr);
      }

      // Add all atoms to cells.
      for (int i = 0; i < nAtom_; ++i) {
         cells_[tags_[i].cellRank].append(tags_[i].handle);
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
         if (handles_.capacity() <= 0) {
            UTIL_THROW("CellList is allocated but handles_.capacity() <= 0");
         }
      }

      if (isBuilt_) {

         if (!isAllocated()) {
            UTIL_THROW("CellList is built but not allocated");
         }

         // Check validity of all cells individually. 
         const Atom* atomPtr;
         const Cell* cellPtr;
         int   nAtomCell;
         int   nAtomSum = 0;
         for (int icell = 0; icell < grid_.size(); ++icell) {
            cellPtr = &cells_[icell];
            nAtomCell = cellPtr->nAtom();
            if (nAtomCell != cellPtr->atomCapacity()) {
               UTIL_THROW("Cell nAtom != atomCapacity");
            }
            if (nAtomCell > 0) {
               for (int i = 0; i < nAtomCell; ++i) {
                  atomPtr = cellPtr->atomPtr(i);
                  if (atomPtr == 0)
                      UTIL_THROW("Null Atom* in a Cell");
                  if (icell != cellIndexFromPosition(atomPtr->position())) {
                      UTIL_THROW("Inconsistent position");
                  }
               }
            }
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
