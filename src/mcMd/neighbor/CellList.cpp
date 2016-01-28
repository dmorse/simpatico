/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CellList.h"
#include <util/space/Vector.h>
#include <util/space/Dimension.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   CellList::CellList()
   {
      for (int i = 0; i < Dimension; ++i) {
         invCellWidths_[i] = 0.0;
         numCells_[i] = 0;
         minCells_[i] = 0;
         maxCells_[i] = 0;
         minDel_[i] = 0;
         maxDel_[i] = 0;
      }
      YZCells_ = 0;
      totCells_ = 0;
      atomCapacity_ = 0;
   }

   /*
   * Destructor. 
   */
   CellList::~CellList()
   {}

   /*
   * Resets all Cell and CellTag objects to empty state.
   */
   void CellList::clear()
   {
      // Clear all Cell objects
      int i;
      if (cells_.capacity() > 0) {
         for (i=0; i < cells_.capacity(); ++i) {
            cells_[i].clear();
         }
      }

      // Clear all CellTag objects
      if (cellTags_.capacity() > 0) {
         for (i = 0; i < cellTags_.capacity(); ++i) {
            cellTags_[i].clear();
         }
      }
   }

   /*
   * Set number of cells and width for one axis (x,y, or z)
   */
   void 
   CellList::setCellsAxis(int axis, double cutoff)
   {
      numCells_[axis] = (int)(lengths_[axis]/cutoff);
      if (numCells_[axis] < 1) {
         numCells_[axis] = 1;
      }
      invCellWidths_[axis] = (double)numCells_[axis];

      minCells_[axis]  = 0;
      maxCells_[axis]  = numCells_[axis]-1;

      if (numCells_[axis] > 2) {
         minDel_[axis] = -1;
         maxDel_[axis] =  1;
      } else if (numCells_[axis] == 2) {
         minDel_[axis] = -1;
         maxDel_[axis] =  0;
      } else if (numCells_[axis] == 1) {
         minDel_[axis] =  0;
         maxDel_[axis] =  0;
      }
   }

   /*
   * Set atomCapacity and, if necessary, allocate cellTags_.
   */
   void 
   CellList::setAtomCapacity(int atomCapacity)
   {
      // Precondition 
      if (atomCapacity <= 0) {
         UTIL_THROW("atomCapacity must be > 0");
      }
      atomCapacity_ = atomCapacity;

      // If necessary, allocate/reallocate cellTags_ array.
      if (cellTags_.capacity() == 0) {
         cellTags_.allocate(atomCapacity_);
      } else 
      if (atomCapacity_ > cellTags_.capacity()) {
         cellTags_.deallocate();
         cellTags_.allocate(atomCapacity_);
      }
   }

   /*
   * Setup an empty cell list.
   */
   void CellList::setup(const Boundary &boundary, double cutoff)
   {
      if (cutoff <= 0) {
         UTIL_THROW("cutoff must be > 0");
      }

      lengths_ = boundary.lengths();
      setCellsAxis(0, cutoff);
      setCellsAxis(1, cutoff);
      setCellsAxis(2, cutoff);
      YZCells_   = numCells_[1]*numCells_[2];
      totCells_  = YZCells_*numCells_[0];
      if (totCells_ < 1) {
         UTIL_THROW("totCells_ must be > 1");
      } 

      // If necessary, allocate or reallocate cells_ array
      if (cells_.capacity() == 0) {
         cells_.allocate(totCells_);
      } else
      if (totCells_ > cells_.capacity()) {
         cells_.deallocate();
         cells_.allocate(totCells_);
      }
      clear();

      boundaryPtr_ = &boundary;
   }

   /*
   * Fill an array with Ids of atoms in cell ic and all neighboring cells.
   */
   void 
   CellList::getCellNeighbors(int ic, NeighborArray &neighbors, int &nInCell)
   const
   {
      const Cell *cellPtr;
      Atom *atomPtr;
      int   icx, icy, icz;
      int   jc, jcx, jcy, jcz;
      int   dcx, dcy, dcz, jp;

      // Determine cell coordinates icx, icy, icz
      cellCoordFromIndex(ic, icx, icy, icz);

      // Zero total number of neighbor atoms and number in cell
      //nNeighbor = 0;
      nInCell   = 0;
      neighbors.clear();

      // Loop over atoms in cell ic.
      // By convention, these appear first in the list.
      cellPtr = &(cells_[ic]);
      for (jp=0; jp < cellPtr->firstClearPos(); ++jp) {
         atomPtr = cellPtr->atomPtr(jp);
         if (atomPtr != 0) {
            neighbors.append(atomPtr);
            ++nInCell;
         }
      }

      // Loop over neighboring cells (excluding cell ic)

      // Loop in x direction
      for (dcx = minDel_[0]; dcx <= maxDel_[0]; ++dcx) {
         jcx = shiftCellCoordAxis(0, icx + dcx);

         // Loop in y direction
         for (dcy = minDel_[1]; dcy <= maxDel_[1]; ++dcy) {
            jcy = shiftCellCoordAxis(1, icy + dcy);

            // Loop in z direction
            for (dcz = minDel_[2]; dcz <= maxDel_[2]; ++dcz) {
               jcz = shiftCellCoordAxis(2, icz + dcz);

               // Get cell index jc and Cell pointer cellPtr
               jc = cellIndexFromCoord(jcx, jcy, jcz);

               if (jc != ic) {

                  // Loop over atoms in neighboring cell jc
                  cellPtr = &cells_[jc];
                  for (jp=0; jp < cellPtr->firstClearPos(); ++jp) {
                     atomPtr  = cellPtr->atomPtr(jp);
                     if (atomPtr != 0) {
                        neighbors.append(atomPtr);
                     }
                  }

               } // end if (jc != ic)

            } // end for dcz

         } // end for dcy

      } // end for dcx

   } 

   /*
   * Get the maximum allowed atom index + 1.
   */
   int CellList::atomCapacity() const
   {  return atomCapacity_; }

   /*
   * Get total number of atoms in this CellList.
   */
   int CellList::nAtom() const
   {
      int nAtomSum = 0;
      for (int icell = 0; icell < totCells_; ++icell) {
         nAtomSum += cells_[icell].nAtomCell();
      }
      return nAtomSum;
   }

   /*
   * Check validity of CellList, throw an Exception if an error is found.
   *
   * This method checks consistency of the CellList and Atom data structures.
   * Calls Cell::isValid internally for each cell. 
   */
   bool CellList::isValid(int nAtom) const
   {

      // Call Cell::isValid for each Cell. Let any Exceptions bubble up.
      int nAtomSum = 0;
      for (int icell = 0; icell < totCells_; ++icell) {
         cells_[icell].isValid(cellTags_, nAtom, icell);
         nAtomSum += cells_[icell].nAtomCell();
      }

      // Check that total number of atoms in all cells equals nAtom.
      if (nAtom >= 0) {
         if (nAtomSum != nAtom) {
            UTIL_THROW("Number of atoms in all cells != nAtom");
         }
      }

      return true;
   }

}
