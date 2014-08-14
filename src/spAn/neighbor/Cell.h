#ifndef SPAN_CELL_H
#define SPAN_CELL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/neighbor/CellAtom.h>
#include <spAn/chemistry/Atom.h>
#include <util/containers/FSArray.h>
#include <util/containers/FArray.h>
#include <util/global.h>

#include <utility>

namespace SpAn
{

   using namespace Util;

   /**
   * A single cell in a CellList.
   *
   * An initialized Cell has an array of CellAtom objects, a pointer to 
   * the next Cell in a linked list, and a pointer to an array of integer
   * offsets to neighboring cells. 
   *
   * A linked list of cells is created by a parent CellList. The method
   * CellList::begin() returns a pointer to the first Cell in the list.
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
   *   Cell::NeighborArray neighbors;    // array of CellAtom* pointers
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
   * \ingroup SpAn_Neighbor_Module
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
      * An array of offsets to the neighboring cells surrounding this cell.
      *
      * Every cell has a pointer to an OffsetArray, which uses relative 
      * cell indices (offsets relative to the cell id of the primary cell) 
      * to identify neighboring cells. There are a maximum of 3 pairs, one
      * for each dimension.  The pair.first contains the offset to move in
      * the positive direction while pair.second contains the offset to move
      * in the negative direction.  The contents of the OffsetArray are
      * calculated in the CellList::makeGrid() function.
      */
      typedef FSArray< std::pair<int,int>, Dimension> OffsetArray;

      /**
      * Static array for holding neighbors in a cell list.
      */
      typedef FSArray<CellAtom*, MaxNeighborAtom> NeighborArray;

      /**
      * Constructor.
      */
      Cell();

      /**
      * Destructor.
      */
      ~Cell();

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
      */
      void setOffsetArray(OffsetArray& offsets);

      /**
      * Reset to empty before incrementing capacity.
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
      * Associate the Cell with an array of CellAtom objects.
      *
      * The final capacity of the cell must be known when this method
      * is called. It associate the Cell with a C array of capacity 
      * Atom* pointers, starting at position begin. It returns a 
      * pointer to an element one past the end of this array segment. 
      *
      * \param begin first element in associated array segment.
      * \return end of array segment (element one past the end)
      */
      CellAtom* initialize(CellAtom* begin);

      /**
      * Append an Atom to an initialized cell.
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
      * Fill an array with pointers to atoms in a cell and neighboring cells.
      *
      * Upon return, the FSArray neighbors contains pointers to all of the
      * atoms in this cell and neighboring cells.  The first nAtom() elements
      * contain pointers to atoms in this cell. 
      *
      * To avoid double counting of pairs, the method only returns atoms from
      * neighboring local cells with a cell id greater than this->id().
      *
      * \param neighbors Array of pointers to neighbor Atoms
      */
      void getNeighbors(NeighborArray& neighbors) const;

   private:

      /// Pointer to first CellAtom in this cell.
      CellAtom* begin_; 

      /// Pointer to neighbor offset array.
      OffsetArray* offsetsPtr_;

      /// Pointer to next local Cell.
      Cell* nextCellPtr_;

      /// Number of atoms in this cell.
      int  nAtom_;

      /// Maximum number of atoms in cell.
      int  atomCapacity_;

      /// Id of cell in grid.
      int id_;

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

}
#endif
