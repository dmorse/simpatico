#ifndef DDMD_CELL_CPP
#define DDMD_CELL_CPP

// THIS CLASS IS UNFINISHED, UNTESTED, AND UNUSED

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Cell.h"

namespace DdMd
{

   using namespace Util;

   Cell::Cell()
    : begin_(0),
      offsetsPtr_(0),
      nextCellPtr_(0),
      nAtom_(0),
      atomCapacity_(0),
      isGhostCell_(true)
   {}

   void Cell::setOffsetArray(Cell::OffsetArray& offsets)
   {  offsetsPtr_ = &offsets; }

   void Cell::setIsGhostCell(bool isGhostCell)
   {  isGhostCell_ = isGhostCell; }

   void Cell::setNextCell(Cell& nextCell)
   {  nextCellPtr_ = &nextCell; }

   void Cell::setLastCell()
   {  nextCellPtr_ = 0; }

   /*
   * Fill an array with pointers to atoms in a cell and neighboring cells.
   *
   * Upon return, the NeighborArray neighbors contains pointers to all of
   * the atoms this cell and neighboring cells.  The first nAtom() elements
   * are the atoms in this cell.
   *
   * \param neighbors array of pointers to neighbor Atoms
   */
   void Cell::getNeighbors(NeighborArray &neighbors) const
   {

      // Preconditions
      assert(offsetsPtr_);
      assert(!isGhostCell_);

      const Cell* cellPtr;
      int   ic, ia, nc, na;

      neighbors.clear();

      nc = offsetsPtr_->size();
      for (int ic = 0; ic < nc; ++ic) {
         cellPtr = this + (*offsetsPtr_)[ic];
         if (cellPtr->id() >= id_ || cellPtr->isGhostCell()) {
            na = cellPtr->nAtom();
            for (ia = 0; ia < na; ++ia) {
               neighbors.append(cellPtr->atomPtr(ia));
            }
         }
      }
   }


}
#endif
