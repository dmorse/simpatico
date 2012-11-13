#ifndef DDMD_CELL_H
#define DDMD_CELL_H

// THIS CLASS IS UNFINISHED, UNTESTED, AND UNUSED

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/chemistry/Atom.h>
#include <util/containers/FSArray.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /**
   * A single Cell in a CellList.
   *
   * An initialized Cell has an array of Atom* pointers, a pointer to 
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
      * Maximum allowed number of neighboring cells. 
      */
      static const int MaxNeighborCell = 27;
      
      /**
      * Maximum possible number of atoms in this an neighboring cells.
      */
      static const int MaxNeighborAtom = 2000;

      /**
      * Static array for holding neighbors in a cell list.
      */
      typedef FSArray<Atom*, MaxNeighborAtom> NeighborArray;

      /**
      * An array of integer offsets to neighboring cells.
      */
      typedef FSArray<int, MaxNeighborCell> OffsetArray;

      /**
      * Identifier type for an Atom.
      */
      typedef Atom* Handle;

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
      * Set the pointer to an array of integer offsets.
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
      * Associate the Cell with an array of Atom* pointers.
      *
      * The final capacity of the cell must be known when this method
      * is called. It associate the Cell with a C array of capacity 
      * Atom* pointers, starting at position begin. It returns a 
      * pointer to an element one past the end of this array segment. 
      *
      * \param begin first element in associated array segment.
      * \return end of array segment (element one past the end)
      */
      Atom** initialize(Atom** begin);

      /**
      * Append an Atom* pointer to an initialized cell.
      */
      void append(Atom* atomPtr);

      // Accessors

      int id() const
      {  return id_; }
   
      void setId(int id) 
      {  id_ = id; }

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
      Atom* atomPtr(int i) const;

      /**
      * Return the number of neighboring cells.
      */
      int nNeighborCell() const;

      /**
      * Is this a ghost cell?
      */
      bool isGhostCell() const;

      /**
      * Return a pointer to neighbor cell i.
      */
      const Cell* neighborCellPtr(int i) const;

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
      Atom**  begin_;         

      /// Pointer to neighbor offset array.
      OffsetArray*  offsetsPtr_;

      /// Pointer to next local Cell.
      Cell*  nextCellPtr_;

      /// Number of atoms in this cell.
      int  nAtom_;

      /// Maximum number of atoms in cell.
      int  atomCapacity_;  

      /// Number of neighboring cells.
      int  nNeighborCell_;

      /// Id of cell in grid.
      int id_;

      /// Is this a ghost cell?
      bool isGhostCell_;

   };

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

   inline Atom** Cell::initialize(Atom** begin)
   {
      assert(begin_ == 0);
      assert(nAtom_  == 0);
      assert(atomCapacity_  >= 0);

      begin_ = begin; 
      return (begin_ + atomCapacity_);
   }

   inline void Cell::append(Atom* atomPtr)
   {
      assert(begin_ != 0);
      assert(nAtom_ < atomCapacity_);
      begin_[nAtom_] = atomPtr;
      ++nAtom_;
   }

   /*
   * Return number of atoms in this cell.
   */
   inline int Cell::nAtom() const
   {  return nAtom_; }

   /*
   * Return pointer to atom i.
   */
   inline Atom* Cell::atomPtr(int i) const
   {
      assert(i < nAtom_);
      return begin_[i];
   }

   /*
   * Return number of neighboring cells. 
   */
   inline int Cell::nNeighborCell() const
   {
      assert(offsetsPtr_);  
      return offsetsPtr_->size(); 
   }

   /*
   * Pointer to neighboring cell number i. 
   */
   inline const Cell* Cell::neighborCellPtr(int i) const
   { 
      assert(offsetsPtr_);  
      return (this + (*offsetsPtr_)[i]); 
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
