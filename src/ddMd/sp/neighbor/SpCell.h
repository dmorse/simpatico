#ifndef DDMD_SP_CELL_H
#define DDMD_SP_CELL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/sp/neighbor/SpCellAtom.h>
#include <ddMd/sp/chemistry/SpAtom.h>
#include <util/containers/FSArray.h>
#include <util/containers/FArray.h>
#include <util/global.h>

#include <utility>

namespace DdMd
{

   using namespace Util;

   /**
   * A single cell in a SpCellList.
   *
   * An initialized SpCell has an array of SpCellAtom objects, a pointer to 
   * the next SpCell in a linked list, and a pointer to an array of integer
   * offsets to neighboring cells. 
   *
   * A linked list of cells is created by a parent SpCellList. The method
   * SpCellList::begin() returns a pointer to the first SpCell in the list.
   *
   * The method SpCell::getNeighbors() returns an array containing pointers
   * to atoms in this cell and all neighboring cells, with the atoms in
   * in this cell listed first.
   *
   * Here is an example of code to iterate over all local cells in a
   * SpCellList, and over all pairs of neighboring atoms:
   * \code
   * 
   *   SpCellList cellList;                // parent SpCellList
   *   SpCell::NeighborArray neighbors;    // array of SpCellAtom* pointers
   *   const SpCell* cellPtr;              // pointer to SpCell in linked list
   *   SpAtom*  atom1Ptr;                  // pointer to Atom in this cell.
   *   SpAtom*  atom2Ptr;                  // pointer to neighbor Atom
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
   class SpCell
   {

   public:

      // Static members

      /**
      * Maximum possible number of atoms in this an neighboring cells.
      */
      static const int MaxNeighborAtom = 2000;

      /**
      * An array of relative cell ids for neighboring cells.
      */
      typedef FArray<int, 27> OffsetArray;

      /**
      * Static array for holding neighbors in a cell list.
      */
      typedef FSArray<SpCellAtom*, MaxNeighborAtom> NeighborArray;

      /**
      * Constructor.
      */
      SpCell();

      // Linked List Interface

      /**
      * Set the pointer to the next cell in the list.
      */
      void setNextCell(SpCell& nextSpCell);

      /**
      * Set this to be the last cell in the list.
      */
      void setLastCell();

      /**
      * Return a pointer to neighbor cell i.
      */
      const SpCell* nextCellPtr() const;

      // Mutators

      /**
      * Set id for this SpCell.
      *
      * \param id integer identifier for this SpCell
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
      * Associate the SpCell with an array of SpCellAtom objects.
      *
      * The final capacity of the cell must be known when this method
      * is called. It associate the SpCell with a C array of capacity 
      * Atom* pointers, starting at position begin. It returns a 
      * pointer to an element one past the end of this array segment. 
      *
      * \param begin first element in associated array segment.
      * \return end of array segment (element one past the end)
      */
      SpCellAtom* initialize(SpCellAtom* begin);

      /**
      * Append an Atom to an initialized cell.
      */
      void append(SpAtom* atomPtr);

      // Accessors

      /**
      * Get identifier for this SpCell.
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
      SpCellAtom* atomPtr(int i) const;

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

      /// Pointer to first SpCellAtom in this cell.
      SpCellAtom* begin_; 

      /// Pointer to neighbor offset array.
      OffsetArray* offsetsPtr_;

      /// Pointer to next local SpCell.
      SpCell* nextCellPtr_;

      /// Number of atoms in this cell.
      int  nAtom_;

      /// Maximum number of atoms in cell.
      int  atomCapacity_;

      /// Id of cell in grid.
      int id_;

   };

   inline void SpCell::setId(int id) 
   {  id_ = id; }

   inline void SpCell::incrementCapacity()
   {
      assert(begin_ == 0);
      ++atomCapacity_;
   }

   inline void SpCell::clear()
   {
      begin_= 0;
      nAtom_ = 0;
      atomCapacity_ = 0;
   }

   inline SpCellAtom* SpCell::initialize(SpCellAtom* begin)
   {
      assert(begin_ == 0);
      assert(nAtom_ == 0);
      assert(atomCapacity_  >= 0);

      begin_ = begin; 
      return (begin_ + atomCapacity_);
   }

   inline void SpCell::append(SpAtom* atomPtr)
   {
      assert(begin_ != 0);
      assert(nAtom_ < atomCapacity_);
      begin_[nAtom_].setPtr(atomPtr);
      ++nAtom_;
   }

   /*
   * Get identifier for this SpCell.
   */
   inline int SpCell::id() const
   {  return id_; }
  
   /*
   * Return number of atoms in this cell.
   */
   inline int SpCell::nAtom() const
   {  return nAtom_; }

   /*
   * Return pointer to atom i.
   */
   inline SpCellAtom* SpCell::atomPtr(int i) const
   {
      assert(i >= 0);
      assert(i < nAtom_);
      return &begin_[i];
   }

   /*
   * Pointer to next cell in list.
   */
   inline const SpCell* SpCell::nextCellPtr() const
   {  return nextCellPtr_; }

   /*
   *  Return current capacity of cell. 
   */
   inline int SpCell::atomCapacity() const
   {  return atomCapacity_; }

}
#endif
