#ifndef DDMD_CELL_H
#define DDMD_CELL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/neighbor/CellAtom.h>
#include <ddMd/chemistry/Atom.h>
#include <util/containers/FSArray.h>
#include <util/global.h>

#include <utility>

namespace DdMd
{

   using namespace Util;

   /**
   * A single Cell in a CellList.
   *
   * An initialized Cell has an array of CellAtom objects, a pointer to 
   * the next Cell in a linked list, and a pointer to an array of integer
   * offsets to neighboring cells. 
   *
   * A linked list of cells is created by a parent CellList. The method
   * CellList::begin() returns a pointer to the first Cell in the list.
   * This linked list normally contains only the cells with local atoms,
   * and excludes cells of ghost atoms.
   *
   * The method Cell::getNeighbors() returns an array containing pointers
   * to atoms in this cell and all neighboring cells, with the atoms in
   * in this cell listed first.
   *
   * Here is an example of code to iterate over all local cells in a
   * CellList, and over all pairs of neighboring atoms:
   * \code
   * 
   *   CellList cellList;                // parent CellList
   *   Cell::NeighborArray neighbors;    // array of Atom* pointers
   *   const Cell* cellPtr;              // pointer to Cell in linked list
   *   Atom*  atom1Ptr;                  // pointer to Atom in this cell.
   *   Atom*  atom2Ptr;                  // pointer to neighbor Atom
   *  
   *   // Iterate over cells in list.
   *   cellPtr = cellList.begin();
   *   while (cellPtr) {
   * 
   *      cellPtr->getNeighbors(neighbors);
   * 
   *      // Iterate over atoms in this cell.
   *      for (i = 0; i < cellPtr->nAtom(); ++i) {
   *         atom1Ptr = neighbors[i];
   * 
   *         // Iterate over neighbor atoms.
   *         for (j = 0; j < neighbors.size(); ++j) {
   *            atom2Ptr = neighbors[j];
   * 
   *         }
   *      }
   *   }
   *
   * \endcode
   * 
   * \ingroup DdMd_Neighbor_Module
   */
   class Cell
   {

   public:

      // Static members

      /**
      * Maximum possible number of atoms in this an neighboring cells.
      */
      static const int MaxNeighborAtom = 2000;

      /**
      * Maximum number of cell per cutoff length.
      */
      static const int MaxNCellCut = 4;
      
      /**
      * Maximum allowed number of neighboring cells. 
      */
      static const int OffSetArrayCapacity = (2*MaxNCellCut + 1)*(2*MaxNCellCut + 1) + 3;
    
      /**
      * An array of strips of relative ids for columns of neighboring cells.
      *
      * Every cell has a pointer to a Cell::OffsetArray, which uses relative 
      * cell indices (offsets relative to the cell id of the primary cell) 
      * to identify neighboring cells. Each std::pair<int, int> element in 
      * in an OffsetArra contains relative addresses for the the first cell
      * (pair.first) and last cell (pair.second) cell in a contiguous strip 
      * of cells that could contain atoms that lie within a cutoff length 
      * of some point in the primary cell. The contents of the OffsetArray
      * are calculated in the CellList::makeGrid() function.
      */
      typedef FSArray< std::pair<int,int>, OffSetArrayCapacity> OffsetArray;

      /**
      * Static array for holding neighbors in a cell list.
      */
      typedef FSArray<CellAtom*, MaxNeighborAtom> NeighborArray;

      /**
      * Constructor.
      */
      Cell();

      // Linked List Interface

      /**
      * Set the pointer to the next cell in the list.
      */
      void setNextCell(Cell& nextCell);

      /**
      * Set this to be the last cell in the list.
      */
      void setLastCell();

      /**
      * Return a pointer to neighbor cell i.
      */
      const Cell* nextCellPtr() const;

      // Mutators 

      /**
      * Set id for this Cell.
      *
      * \param id integer identifier for this Cell
      */
      void setId(int id);

      /**
      * Set the pointer to an array of integer offsets.
      *
      * \param offsets array offsets for strips of neighboring cells
      */
      void setOffsetArray(OffsetArray& offsets);

      /**
      * Mark as a ghost or local cell.
      */
      void setIsGhostCell(bool isGhostCell = true);

      /**
      * Reset to empty before incrementing capacity.
      *
      * Does not nullify nextCell pointer is isGhostCell bool.
      */
      void clear();

      /**
      * Increment the capacity counter.
      *
      * This must be called within a loop over atoms, once
      * per atom that belongs in this cell. This loop must
      * completed before initialize is called.
      */
      void incrementCapacity();

      /**
      * Associate this Cell with an array block of CellAtom objects.
      *
      * The final atomCapacity of the cell must be known when this function
      * is called. It associates the Cell with a C array of capacity Atom*
      * pointers, with the first element at position given by the pointer
      * parameter "begin". This function returns a pointer to an element 
      * one past the end of the associated array segment, which can be used 
      * as the first element of the next Cell. 
      *
      * \param begin first element in associated array segment.
      * \return end of array segment (element one past the end)
      */
      CellAtom* initialize(CellAtom* begin);

      /**
      * Append an Atom to an initialized cell.
      *
      * This function should be called once per atom, after the cell has
      * been initialized by calling CellAtom::initialize.
      */
      void append(Atom* atomPtr);

      // Accessors

      /**
      * Get identifier for this Cell.
      */
      int id() const;
  
      /**
      * Number of atoms in cell.
      */
      int nAtom() const;

      /**
      * Capacity of array segment. 
      */
      int atomCapacity() const;

      /**
      * Return a pointer to atom i.
      */
      CellAtom* atomPtr(int i) const;

      /**
      * Is this a ghost cell?
      */
      bool isGhostCell() const;

      /**
      * Fill an array with pointers to atoms in a cell and neighboring cells.
      *
      * Upon return, the FSArray neighbors contains pointers to all of the
      * atoms in this cell and neighboring cells.  The first nAtom() elements
      * contain pointers to atoms in this cell. 
      *
      * To avoid double counting of pairs, the method only returns atoms from
      * neighboring local cells with a cell id greater than this->id(), and 
      * from neighboring ghost cells.
      *
      * \param neighbors          Array of pointers to neighbor Atoms
      * \param reverseUpdateFlag  Is reverse communication enabled?
      */
      void getNeighbors(NeighborArray& neighbors, 
                        bool reverseUpdateFlag = false) const;

   private:

      /// Pointer to first Atom* pointer for this cell.
      CellAtom* begin_;         

      /// Pointer to neighbor offset array.
      OffsetArray* offsetsPtr_;

      /// Pointer to next local Cell.
      Cell* nextCellPtr_;

      /// Number of atoms currently in this cell.
      int nAtom_;

      /// Number of atoms for which space is assigned.
      int atomCapacity_;  

      /// Id of cell in grid.
      int id_;

      /// Is this a ghost cell?
      bool isGhostCell_;

      // Note: The atomCapacity_ should be incremented during the initial
      // loop over atoms by calling Cell::incrementCapacity() each time 
      // an atom is found to be in this Cell. The nAtom_ member is 
      // incremented by the Cell:append method each time an atoms is
      // actually added to this cell.

   };

   inline void Cell::setId(int id) 
   {  id_ = id; }

   inline void Cell::incrementCapacity()
   {
      assert(begin_ == 0);
      ++atomCapacity_;
   }

   inline void Cell::clear()
   {
      begin_= 0;
      nAtom_ = 0;
      atomCapacity_ = 0;
   }

   inline CellAtom* Cell::initialize(CellAtom* begin)
   {
      assert(begin_ == 0);
      assert(nAtom_ == 0);
      assert(atomCapacity_  >= 0);

      begin_ = begin; 
      return (begin_ + atomCapacity_);
   }

   inline void Cell::append(Atom* atomPtr)
   {
      assert(begin_ != 0);
      assert(nAtom_ < atomCapacity_);
      begin_[nAtom_].setPtr(atomPtr);
      ++nAtom_;
   }

   /*
   * Get identifier for this Cell.
   */
   inline int Cell::id() const
   {  return id_; }
  
   /*
   * Return number of atoms in this cell.
   */
   inline int Cell::nAtom() const
   {  return nAtom_; }

   /*
   * Return pointer to atom i.
   */
   inline CellAtom* Cell::atomPtr(int i) const
   {
      assert(i >= 0);
      assert(i < nAtom_);
      return &begin_[i];
   }

   /*
   * Pointer to next cell in list.
   */
   inline const Cell* Cell::nextCellPtr() const
   {  return nextCellPtr_; }

   /*
   *  Return current capacity of cell. 
   */
   inline int Cell::atomCapacity() const
   {  return atomCapacity_; }

   /*
   * Is this a ghost cell?
   */
   inline bool Cell::isGhostCell() const
   {  return isGhostCell_; }

}
#endif
