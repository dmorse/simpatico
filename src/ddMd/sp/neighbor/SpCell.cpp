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
   * Upon return, the NeighborArray neighbors contains pointers to all of the
   * atoms this cell and neighboring cells.  The first nAtom() elements are the
   * the atoms in this cell.
   */
   void SpCell::getNeighbors(NeighborArray &neighbors) const
   {
      // Preconditions
      assert(offsetsPtr_);
      assert(!isGhostCell_);

      const SpCell* cellBegin;
      const SpCell* cellEnd;
      SpCellAtom* atomBegin;
      SpCellAtom* atomEnd;
      int  is, ns;
      bool bg, eg;

      neighbors.clear();
      ns = offsetsPtr_->size();

      for (is = 0; is < ns; ++is) {
         cellBegin = this + (*offsetsPtr_)[is].first;
         cellEnd   = this + (*offsetsPtr_)[is].second;
         if (cellBegin->id() >= id_) {
            atomBegin = cellBegin->begin_;
            atomEnd = cellEnd->begin_ + cellEnd->nAtom_;
            for ( ; atomBegin < atomEnd; ++atomBegin) {
               neighbors.append(atomBegin);
            }
         }
      }
   }

}
#endif
