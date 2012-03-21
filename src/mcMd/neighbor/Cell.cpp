#ifndef  CELL_CPP
#define MCMD_CELL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Cell.h"
#include <mcMd/chemistry/Atom.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Default constructor, creates an empty Cell.
   */
   Cell::Cell() 
   {  clear();}
   
   /* 
   * Reset Cell to empty state.
   */
   void Cell::clear() 
   {
      nAtomCell_     = 0;
      firstEmptyPos_ = 0;
      firstClearPos_ = 0;
      for (int j = 0; j < MaxAtomCell; ++j) {
         atoms_[j] = 0;
      }
   }
   
   
   /* 
   * Check validity of Cell data members.
   */
   bool Cell::isValid(const Array<CellTag> &cellTags, int nAtom, int icell) 
   const
   {
      int jp, atomId, nAtomCellTest;
      nAtomCellTest=0;
      for (jp=0; jp < Cell::MaxAtomCell; ++jp) {
         if (atoms_[jp] != 0) {
            atomId = atoms_[jp]->id();
            ++nAtomCellTest;
            if (cellTags[atomId].cellId != icell) {
               UTIL_THROW("Value in CellTag.cellId inconsistent with icell");
            }
            if (cellTags[atomId].cellPos != jp) {
	       std::cout << "atomId                   =" << atomId << std::endl;
	       std::cout << "cellTags[atomId].cellPos =" 
		         << cellTags[atomId].cellPos << std::endl;
	       std::cout << "index jp in Cell         = " <<  jp << std::endl;
               UTIL_THROW(
                     "Value in CellTag.cellPos inconsistent with Cell.atoms_");
            }
            if (jp >= firstClearPos_) {
               UTIL_THROW("Value atoms_[i] >=0 for i >= firstClearPos_");
            }
         } else {
            if (jp < firstEmptyPos_) {
               UTIL_THROW("Value atoms_[i] < 0 for i < firstEmptyPos_");
            }
         }
      }
      if (firstClearPos_ > 0 && atoms_[firstClearPos_-1] < 0) {
         UTIL_THROW("Error: firstClearPos_: atoms_[firstClearPos_-1] < 0");
      }
      if (nAtomCellTest != nAtomCell_) {
         UTIL_THROW("Number of atoms in cell != nAtomCell_");
      }
      return true;
   }

} 
#endif
