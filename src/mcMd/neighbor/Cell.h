#ifndef MCMD_CELL_H
#define MCMD_CELL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CellTag.h" 
#include <util/containers/Array.h> 
#include <util/global.h>

class CellTest;
class CellListTest;

namespace McMd
{

   using namespace Util;

   class Atom;
   
   /**
   * A set of Atoms in a small region.
   *
   * A Cell object holds an array of pointers to the Atoms that lie within 
   * a small region within a System. A CellList contains a private array of
   * Cell objects.  A Cell should be accessed only by its parent CellList. 
   *
   * The maximum allowed number of atoms in a cell is given by a static 
   * constant MaxAtomCell. Attempting to add more than MaxAtomCell atoms 
   * will cause an Exception to be thrown by Cell::addAtom(), thus ending
   * the simulation. If this happens, the only solution is to recompile 
   * with a larger value of MaxAtomCell.
   *
   * \ingroup McMd_Neighbor_Module
   */
   class Cell 
   {
     
   public:

      /// Maximum number of atoms per cell.
      #ifdef SIMP_COULOMB
      static const int MaxAtomCell = 253; 
      #else
      static const int MaxAtomCell = 61; 
      #endif
   
      /// A null (uninitialized) index value.
      static const int NullIndex   = -1; 
   
      /** 
      * Default constructor, creates an empty Cell.
      */
      Cell();
   
      /**
      * Reset Cell to empty state.
      */
      void clear();
   
      /**
      * Remove a specified atom from the Cell.
      *
      * \param cellTag cellTag object associated with atom to be removed.
      */
      void deleteAtom(CellTag &cellTag);
    
      /**
      * Add one atom to a Cell.
      *
      * On return: cellTag.cellPos is set, and cellTag.cellId = cellId.
      *
      * \param cellTag    CellTag objects associated with added atom
      * \param atom       atom to be added
      * \param cellId     integer Id of this Cell object
      */
      void addAtom(CellTag &cellTag, Atom &atom, int cellId);
 
      /// Get number of atoms in cell. 
      int nAtomCell() const;

      /// Get index of first clear element of atoms_ 
      int firstClearPos() const; 
        
      /** 
      * Get a pointer to Atom (may be null). 
      *
      * \param index element index of atoms_ array
      */
      Atom* atomPtr(int index) const;
        
      /**
      * Check validity of Cell, and consistency with array of CellTags.
      *
      * The parameter cellTags is an array of CellTag objects, in which the
      * array index of each CellTag is the Id of the associated atom.
      *
      * This method throws an Exception if an error is found. If it completes
      * normally, than the Cell passed all tests.
      *
      * \param cellTags  array of CellTag objects
      * \param nAtom     number of atoms in array cellTags
      * \param icell     index of this Cell in parent CellList
      * \return true if valid, throw exception otherwise
      */
      bool isValid(const Array<CellTag> &cellTags, int nAtom, int icell) const;
   
   private:
   
      /*
      * Implementation notes:
      *
      * Pointers to atoms in a cell are stored in the atoms_[] array member.
      * "Empty" elements of atoms_[] must contain a null (i.e., 0) pointer.
      *
      * The total number of atoms in the cell in is nAtomCell_. The array 
      * index of the first empty element in the atoms_[] array is given by 
      * firstEmptyPos_. All elements of atoms_[] with indices greater than 
      * or equal to firstClearPos_ are empty. An element atoms_[i] with 
      * i < firstEmptyPos_ may not be empty, one with i >= firstClearPos_ must 
      * be empty, and one with firstEmptyPos_ < i and i < firstClearPos_ may or 
      * may not be empty. When firstClearPos_ > 0, atoms_[firstClearPos_-1] may 
      * not be empty, since this would imply that firstClearPos_ is not the 
      * first in the contiguous block of empty elements. These conditions are 
      * checked in the isValid() method. 
      */
   
      /// C array containing pointers to atoms in this Cell.
      Atom* atoms_[MaxAtomCell]; ///< stores pointers to the atoms in cell
   
      /// Number of atoms currently in this Cell.
      int   nAtomCell_;   
   
      /// Index of the first empty element in array atoms_[].
      int   firstEmptyPos_;  
   
      /// Index of the first element in a contiguous block of empty elements.
      int   firstClearPos_; 

   //friends:

      //  Grant access to unit test classes
      friend class ::CellTest;
      friend class ::CellListTest;

   }; 
  
   // Inline functions
   
   /* 
   * Delete atom from Cell
   */
   inline void Cell::deleteAtom(CellTag &cellTag) 
   {
      int cellPos = cellTag.cellPos;
   
      // Delete atom from Cell member variables.
      // Set atoms_[cellPos] to 0 to mark empty.
      atoms_[cellPos] = 0;
      nAtomCell_--;
   
      // Reset firstEmptyPos_
      if (firstEmptyPos_ > cellPos) {
         firstEmptyPos_ = cellPos;
      }
   
      // Reset firstClearPos_
      int i = firstClearPos_ - 1;
      while (i >= 0 && atoms_[i] == 0) {
         i--;
      }
      firstClearPos_ = i + 1;
   
      // Reset cellTag cellId and cellPos to negative integer values 
      // to indicate that the atom is no longer associated with a Cell.
      cellTag.cellId  = NullIndex;
      cellTag.cellPos = NullIndex;
   
   }
   
   /* 
   * Add a atom to a Cell
   */
   inline void Cell::addAtom(CellTag &cellTag, Atom &atom, int cellId) 
   {
   
      // If cell is already full, throw Exception
      if (firstEmptyPos_ == MaxAtomCell) {
         UTIL_THROW("Too many atoms in one cell");
      }
       
      // Add atom to Cell
      atoms_[firstEmptyPos_] = &atom;
      ++nAtomCell_;
   
      // Record cell in cellTag
      cellTag.cellId  = cellId;
      cellTag.cellPos = firstEmptyPos_;
   
      // Reset firstEmptyPos_ and firstClearPos_
      if (firstEmptyPos_ < MaxAtomCell - 1) {
   
         if (firstEmptyPos_ == firstClearPos_) {
   
            firstEmptyPos_ = firstEmptyPos_ + 1;
            firstClearPos_ = firstClearPos_ + 1;
   
         } else {
   
            // Find the next empty element, for which atoms_[i] is null
            for (int i = firstEmptyPos_ + 1; i < MaxAtomCell; ++i) {
               if (atoms_[i] == 0) {
                  firstEmptyPos_ = i;
                  break;
               }
            }
   
         }
   
      } else { 
   
         // If firstEmptyPos_ was MaxAtomCell - 1, no empty slots remain.
         // Set firstEmptyPos_ = MaxAtomCell to indicate this condition.
         // Attempting to add one more atom will trigger an Exception.
         firstEmptyPos_ = MaxAtomCell;  // Note: Invalid array index
         firstClearPos_ = MaxAtomCell;  // Note: Invalid array index
   
      }
   }

   inline int Cell::nAtomCell() const
   { return nAtomCell_; }
         
   inline int Cell::firstClearPos() const
   { return firstClearPos_; }

   inline Atom* Cell::atomPtr(int index) const
   { return atoms_[index]; }
         
} 

#endif
