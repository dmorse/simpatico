#ifndef DDMD_SP_CELL_CPP
#define DDMD_SP_CELL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SpCell.h"

namespace DdMd
{

   using namespace Util;

   SpCell::SpCell()
    : begin_(0),
      offsetsPtr_(0),
      nextCellPtr_(0),
      nAtom_(0),
      atomCapacity_(0)
   {}

   void SpCell::setOffsetArray(SpCell::OffsetArray& offsets)
   {  offsetsPtr_ = &offsets; }

   void SpCell::setNextCell(SpCell& nextCell)
   {  nextCellPtr_ = &nextCell; }

   void SpCell::setLastCell()
   {  nextCellPtr_ = 0; }

   /*
   * Fill an array with pointers to SpCellAtom objects in this cell and neighbors.
   *
   * Upon return, the NeighborArray neighbors contains pointers to all of the atoms
   * in this cell and neighboring cells.  The first nAtom() elements are the the 
   * atoms in this cell.
   */
   void SpCell::getNeighbors(NeighborArray &neighbors) const
   {
      // Preconditions
      assert(offsetsPtr_);

      const SpCell* cellPtr;
      SpCellAtom* atomBegin;
      SpCellAtom* atomEnd;

      neighbors.clear();

      for (int is = 0; is < 27; ++is) {
         cellPtr = this + (*offsetsPtr_)[is];
         if (cellPtr->id() >= id_) {
            atomBegin = cellPtr->begin_;
            atomEnd = atomBegin + cellPtr->nAtom_;
            for ( ; atomBegin < atomEnd; ++atomBegin) {
               neighbors.append(atomBegin);
            }
         }
      }
   }

}
#endif
