#ifndef MCMD_CELL_LIST_H
#define MCMD_CELL_LIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Cell.h"
#include "CellTag.h"
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Atom.h>
#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/global.h>

#include <sstream>

class CellListTest;

namespace McMd
{

   using namespace Util;

   /**
   * A cell list for Atom objects in a periodic system boundary. 
   *
   * A CellList divides the volume within a Boundary into a grid of 
   * cells (i.e., regions), such that the length of each cell in each
   * of 3 directions is greater * than a specified cutoff distance. 
   * Each cell has an integer index between 0 and totCells_ - 1. 
   * A CellList can be used to identify all atoms that are within a 
   * cutoff distance of any point within a system, as required to 
   * calculate nonbonded energies and forces.
   *
   * In this implementation, each cell is represented by a Cell
   * object. Each cell contains an small array of pointers to all
   * all atoms in the associated volume. This implementation is
   * designed to allow fast addition and removal of individual 
   * atoms. This implementation also wastes some memory, compared
   * to a linked list implementation, because the dimension
   * Cell::MaxAtomCell of the array of pointers in each Cell must
   * be large enough to accomodate the maximum number of atoms 
   * that will ever be encountered in an individual cell.
   *
   * CellList is a non-polymorphic class, with no virtual methods, 
   * and a non-virtual destructor. Do not derive subclasses from it.
   *
   * \ingroup McMd_Neighbor_Module
   */
   class CellList
   {

   public:

      // Static members

      /**
      * Maximum possible number of neighboring atoms.
      *
      * Use as dimension of the array neighbor passed to getNeighbors() and
      * getCellNeighbors().
      */
      static const int MaxNeighbor = 27*Cell::MaxAtomCell;

      /**
      * Static array for holding neighbors in a cell list.
      *
      * \ingroup McMd_Neighbor_Module
      */
      typedef FSArray<Atom*, MaxNeighbor> NeighborArray;

      // public methods

      /**
      * Constructor.
      */
      CellList();

      /**
      * Destructor.
      *
      * Deallocates array of Cell objects, if necessary.
      */
      virtual ~CellList();

      /**
      * Allocate memory for this CellList.
      *
      * This function:
      *   - Allocates an array of atomCapacity CellTag objects.
      *   - Allocates an array of Cell objects sized for this boundary.
      *
      * The capacity of the array of Cells is chosen such that, if
      * each dimension of the system increases by the parameter 
      * "scale", all dimensions of each cell will remain greater than the 
      * parameter "cutoff". The cutoff should be greater than or equal to 
      * the maximum range of the nonbonded pair potential. 
      *
      * This function does not populate the CellList with atoms. See the
      * CellList::addAtom() methods.
      *
      * \param atomCapacity dimension of global array of atoms
      * \param boundary     maximum Boundary for allocation of array of cells
      * \param cutoff       minimum dimension of a cell in any direction
      */
      void allocate(int atomCapacity, const Boundary &boundary, double cutoff);

      /**
      * Initialize grid geometry.
      *
      * The number of cells in each direction is chosen such that the dimension
      * of each cell in each direction is greater than or equal to the cutoff
      * parameter. To calculate nonbonded pair interaction energies, the cutoff
      * parameter should thus be equal to or greater than the maximum range of
      * nonbonded interactions.
      *
      * \param boundary Boundary object for the system
      * \param cutoff   Minimum dimension of a cell in any direction
      */
      void makeGrid (const Boundary &boundary, double cutoff);

      /**
      * Sets all Cell objects to empty state (no Atoms).
      */
      void clear();

      /**
      * Return index of the cell that contains the position pos = (x,y,z).
      *
      * \param pos position vector, as an array pos = {x, y, z}
      * \return index of cell that contains position pos[3]
      */
      int cellIndexFromPosition(const Vector &pos) const;

      /**
      * Add a Atom to the appropriate cell, based on its position.
      *
      * \param atom  Atom object to be added.
      */
      void addAtom(Atom &atom);

      /**
      * Delete a Atom object from its cell.
      *
      * \param atom Atom object to be deleted
      */
      void deleteAtom(Atom &atom);

      /**
      * Update CellList to reflect a new atom position.
      *
      * If a new atom position pos[3] lies in a new cell, with a cell index
      * icell != atom.cellId, then delete the atom from the old Cell
      * object and add it to the new one. If the new cell is the same as the
      * old one (icell == atom.CellId), do nothing and return.
      *
      * This method updates the CellList, but does not update the coordinates
      * stored in the pos member of the associated Atom object. To update
      * both, use the System::moveAtom() method.
      *
      * \param atom Atom object to be moved
      * \param pos  new Atom position
      */
      void updateAtomCell(Atom &atom, const Vector &pos);

      /**
      * Fill a NeighborArray with pointers to atoms near a specified position.
      *
      * Upon return, array neighbors contains pointers to all Atom objects 
      * in the cell containing the position pos, and in all neighboring cells.
      * The number of such Atoms is given by the size() of this FSArray < Atom* >
      * array.
      *
      * \param pos        position Vector, pos = {x, y, z}
      * \param neighbors  array of pointers to neighbor Atoms.
      */
      void
      getNeighbors(const Vector &pos, NeighborArray &neighbors) const;

      /**
      * Fill an array with pointers to atoms in a cell and neighboring cells.
      *
      * Upon return, the NeighborArray neighbors contains pointers to all of
      * the atoms in cell number ic and all neighboring cells.  The number 
      * of such Atoms is given by the size() of this FSArray < Atom* > array.
      * Also, upon return, nInCell is the number of Atoms in cell ic. 
      * Pointers to the atoms in cell ic are listed in the elements
      * neighbors[0], ... , neighbors[nInCell-1] of the neighbors array, 
      * while those in neighboring are in elements neighbors[nInCell], 
      * ... , neighbors[neighbors.size()-1].
      *
      * \param ic        cell index
      * \param neighbors FSArray containing pointers to neighbor Atoms
      * \param nInCell   number of Atoms in cell ic.
      */
      void 
      getCellNeighbors(int ic, NeighborArray &neighbors, int &nInCell) const;

      /**
      * Number of cells along axis i.
      *
      * \param i index for axis (direction).
      */
      int gridDimension(int i) const;

      /**
      * Get total number of cells in this CellList.
      */
      int totCells() const;

      /**
      * Get total number of atoms in this CellList.
      */
      int nAtom() const;

      /**
      * Has memory been allocated for this CellList?
      */
      bool isAllocated() const;

      /**
      * Return true if valid, or throw Exception.
      *
      * If nAtom >= 0, check that number of atoms equals nAtom.
      * If nAtom < 0, skip this check.
      *
      * \param nAtom  number of atoms in atom and position arrays.
      * \return true if valid, otherwise throw an Exception.
      */
      bool isValid(int nAtom=-1) const;

   private:

      /// Array of Cell objects
      DArray<Cell>    cells_;

      /// Array of CellTag objects for quick retrieval
      DArray<CellTag> cellTags_;

      Vector lengths_;         ///<  boundary lengths in each direction
      Vector invCellWidths_;   ///<  dimensions of each Cell in CellList grid
      int    minDel_[3];       ///<  Min difference along each axis (usually -1)
      int    maxDel_[3];       ///<  Max difference along each axis (usually +1)
      int    numCells_[3];     ///<  Number of cells in each direction
      int    minCells_[3];     ///<  Min index of cell for each axis = 0
      int    maxCells_[3];     ///<  Max cell index for each axis = numCells - 1
      int    YZCells_;         ///<  numCells_[1]*numCells_[2]
      int    totCells_;        ///<  total number of cells in CellList

      /// Pointer to associated Boundary (set in makeGrid).
      const Boundary* boundaryPtr_; 

      /**
      * Returns the number of cells required along a single Cartesian axis.
      *
      * The number of cells is chosen so that every dimension of each cell 
      * is greater than or equal to parameter cutoff. The pointer member
      * boundary must point to a valid Boundary object.
      *
      * \param boundary  Boundary object containing system dimensions
      * \param axis      index of a Cartesian direction. Values: 0, 1, or 2
      * \param cutoff    minimum allowed dimension for a cell
      */
      void setCellsAxis(const Boundary &boundary, int axis, double cutoff);

      /**
      * Return shifted integer cell coordinate x for axis i.
      *
      * Return integer cell coordinate x for axis i shifted to the range
      * minCells_[i] <= x <= maxCells_[i].
      *
      * \param i index for Cartesian axis. Must be 0, 1, or 2.
      * \param x integer coordinate of cell along axis i.
      */
      int shiftCellCoordAxis(int i, int x) const;

      /**
      * Return integer index of a cell with integer coordinates (cx,cy,cz).
      *
      * \param cx  coordinate of cell along x direction.
      * \param cy  coordinate of cell along y direction.
      * \param cz  coordinate of cell along z direction.
      */
      int cellIndexFromCoord(int cx, int cy, int cz) const;

      /**
      * Return (by reference) int coordinates (cx,cy,cz) of cell with index i.
      */
      void cellCoordFromIndex(int i, int &cx, int &cy, int &cz) const;

      /**
      * Return true if atomId is valid, i.e., if 0 <= 0 < atomCapacity.
      */
      bool isValidAtomId(int atomId)
      { return ( (0 <= atomId) && (atomId < cellTags_.capacity()) ); }


   //friends:

      /// Grant access to unit test class.
      friend class ::CellListTest;

   }; // end class CellList


   // Public inline method definitions:

   /*
   * Get total number of cells in this CellList.
   */
   inline int CellList::totCells() const
   { return totCells_; }

   /*
   *  Return index of the cell that contains the position array pos = {x, y, z}.
   */
   inline int CellList::cellIndexFromPosition(const Vector &pos) const
   {
      int cx, cy, cz;

      if (UTIL_ORTHOGONAL) {
         cx = int(pos[0]*invCellWidths_[0]);
         cy = int(pos[1]*invCellWidths_[1]);
         cz = int(pos[2]*invCellWidths_[2]);
      } else {
         Vector posG;
         boundaryPtr_->transformCartToGen(pos, posG);
         cx = int(posG[0]*invCellWidths_[0]);
         cy = int(posG[1]*invCellWidths_[1]);
         cz = int(posG[2]*invCellWidths_[2]);
      }

      #if 0
      // need to handle case where wrong particle bin is assigned
      // to particles close to the boundary due to roundoff error
      if (cx == numCells_[0]) {
         if (pos[0] >= lengths_[0])
            UTIL_THROW("Particle left boundary.");
         cx = numCells_[0] - 1;
      }
      if (cy == numCells_[1]) {
         if (pos[1] >= lengths_[1])
            UTIL_THROW("Particle left boundary.");
         cy = numCells_[1] - 1;
      }
      if (cz == numCells_[2]) {
         if (pos[2] >= lengths_[2])
            UTIL_THROW("Particle left boundary.");
         cz = numCells_[2] - 1;
      }
      #endif

      assert(cx >= 0);
      assert(cx < numCells_[0]);
      assert(cy >= 0);
      assert(cy < numCells_[1]);
      assert(cz >= 0);
      assert(cz < numCells_[2]);

      return cz + cy*numCells_[2] + cx*YZCells_;
   }

   /*
   * Remove an Atom object from its Cell.
   */
   inline void CellList::deleteAtom(Atom &atom)
   {
      int atomId   = atom.id();
      assert(isValidAtomId(atomId));
      CellTag& cellTag = cellTags_[atomId];
      cells_[cellTag.cellId].deleteAtom(cellTag);
   }

   /*
   * Add a Atom to the appropriate cell, based on its position.
   */
   inline void CellList::addAtom(Atom &atom)
   {
      Vector position = atom.position();
      int    cellId   = cellIndexFromPosition(position);
      int    atomId   = atom.id();
      assert(isValidAtomId(atomId));
      cells_[cellId].addAtom(cellTags_[atomId], atom, cellId);
   }

   /*
   * Update CellList to reflect new atom position
   */
   inline void CellList::updateAtomCell(Atom &atom, const Vector &pos)
   {
      int atomId = atom.id();
      assert(isValidAtomId(atomId));
      CellTag& cellTag = cellTags_[atomId];
      int      oldCell = cellTag.cellId;
      int      newCell = cellIndexFromPosition(pos);
      if (oldCell != newCell) {
         cells_[oldCell].deleteAtom(cellTag);
         cells_[newCell].addAtom(cellTag, atom, newCell);
      }
   }

   /*
   * Fill an array with Ids of atoms that are near a specified position.
   */
   inline void
   CellList::getNeighbors(const Vector &pos, NeighborArray &neighbors) const
   {
      int nInCell;
      int ic = cellIndexFromPosition(pos);
      getCellNeighbors(ic, neighbors, nInCell);
   }

   // Private inline methods:

   /*
   * Return value of integer cell coordinate x for axis i shifted to the range
   * minCells_[i] <= x <= maxCells_[i].
   */
   inline int CellList::shiftCellCoordAxis(int i, int x) const
   {
      if (x > maxCells_[i]) return x - numCells_[i];
      if (x < minCells_[i]) return x + numCells_[i];
      return x;
   }

   /*
   * Return integer index of a cell with integer coordinates (cx,cy,cz).
   */
   inline int CellList::cellIndexFromCoord(int cx, int cy, int cz) const
   {
      return cx*YZCells_ + cy*numCells_[2] + cz ;
   }

   /*
   * Return (by reference) int coordinates (cx,cy,cz) of cell with index i.
   */
   inline
   void CellList::cellCoordFromIndex(int i, int &cx, int &cy, int &cz) const
   {
      cx = (int)(i / YZCells_);
      i = i - cx*YZCells_;
      cy = (int)(i / numCells_[2]);
      cz = i - cy*numCells_[2];
   }

   inline int CellList::gridDimension(int i) const
   {  return numCells_[i]; }

   inline bool CellList::isAllocated() const
   {  return (cells_.capacity() > 0); }

} 
#endif
