#ifndef DDMD_CELL_LIST_H
#define DDMD_CELL_LIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Cell.h"
#include <ddMd/chemistry/Atom.h>
#include <simp/boundary/Boundary.h>
#include <util/space/Grid.h>
#include <util/containers/DArray.h>
#include <util/containers/GArray.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * A cell list used only to identify nearby atom pairs.
   *
   * An CellList divides the domain owned by a processor, plus a frame
   * containing ghost particles, into a grid of cells, such that the length
   * of each cell in each direction is greater than a specified cutoff
   * distance. The algorithm works with either generalized or Cartesian
   * coordinates, if used consistently. The usage example given below
   * assumes atomic coordinates are in generalized/scaled form on entry.
   *
   * All operations of this class are local (no MPI).
   *
   * Usage: 
   * \code
   *
   *    AtomStorage storage;
   *    Boundary boundary;
   *    CellList cellList;
   *    Vector   lower;        // Vector of lower bounds (local atoms)
   *    Vector   upper;        // Vector of upper bounds (local atoms)
   *    Vector   cutoffs;      // Vector of cutoff lengths for each axis.
   *    int      atomCapacity  // max number of atoms on this processor
   *
   *    // Allocate space for structures that store atom info
   *    cellList.setAtomCapacity(atomCapacity);
   *  
   *    // Set elements of cutoffs vector (scaled coordinate)
   *    for (int i = 0; i < Dimension; ++i) {
   *       cutoffs[i] = cutoff/boundary.length(i);
   *    } 
   *
   *    // Setup the grid of cells, and clear all cells
   *    cellList.makeGrid(lower, upper, cutoffs);
   *    cellList.clear();
   *
   *    // Place all local atoms (compute and store cell indices)
   *    AtomStorage::AtomIterator  atomIter;
   *    for (storage.begin(atomIter); atomIter.notEnd(); ++atomIter) {
   *       cellList.placeAtom(*atomIter);
   *    }
   *
   *    // Place all ghost atoms
   *    AtomStorage::GhostIterator ghostIter;
   *    for (storage.begin(ghostIter); ghostIter.notEnd(); ++ghostIter){
   *       cellList.placeAtom(*ghostIter);
   *    }
   *    
   *    // Build cell list
   *    cellList.build();
   *
   *    // Transform atomic positions back to Cartesian coordinates
   *    storage.transformGenToCart(boundary);
   *
   *    // Update atomic positions stored in CellList.
   *    cellList.update();
   *
   * \endcode
   *
   * The atomCapacity parameter should be set equal to the sum of atomCapacity 
   * and the ghostCapacity of the associated atomStorage, which is the maximum 
   * total number of atoms that can exist on this processor.
   *
   * If the upper and lower bounds and atom coordinates are all expressed in
   * generalized coordinates, which span [0,1] over the primitive periodic
   * cell in each direction, each element of the cutoffs vector is given by
   * a ratio cutoffs[i] = cutoff/length[i], where length[i] is the Cartesian
   * distance across the unit cell along a direction parallel to reciprocal 
   * basis vector i, or perpendicular to a surface of the unit cell.
   *
   * See Cell documentation for an example of how to iterate over local cells 
   * and neighboring atom pairs. 
   * 
   * \ingroup DdMd_Neighbor_Module
   */
   class CellList
   {

   public:

      // Public methods

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
      * Set atomCapacity, and allocate arrays that hold atoms.
      *
      * This function:
      *
      *   - Allocates an array of atomCapacity CellList::Tag objects.
      *   - Allocates an array of atomCapacity CellAtom objects.
      *
      * \param atomCapacity maximum number of atoms on this processor
      */
      void setAtomCapacity(int atomCapacity);

      /**
      * Make the cell grid (using generalized coordinates).
      *
      * This method makes a Cell grid in which the number of cells in each
      * direction i is chosen such that the dimension of each cell that 
      * direction is greater than or equal to cutoff[i].
      *
      * The elements of lower and upper should be upper and lower bounds 
      * for coordinates of local atoms on this processor, in generalized
      * coordinates in which the entire periodic unit cell spans 0.0 to 1.0.
      * On input, each element of cutoff[i] should be equal to the minimum
      * cell length in direction i in generalized coordinates. This is given
      * by the ratio cutoff[i] = pairCutoff/length[i], where pairCutoff is 
      * the maximum range of nonbonded interactions, and length[i] is the 
      * distance across the primitive unit cell along the direction 
      * parallel to reciprocal lattice basis vector i.
      *
      * \param lower  lower bound of local atom coordinates.
      * \param upper  upper bound of local atom coordinates.
      * \param cutoffs  pair cutoff length in each direction
      * \param nCellCut  number of cells per cutoff length
      */
      void makeGrid(const Vector& lower, 
                    const Vector& upper, const Vector& cutoffs, 
                    int nCellCut = 1);

      /**
      * Determine the appropriate cell for an Atom, based on its position.
      *
      * This function should be called within a loop over all atoms in the
      * system, which should be completed before invoking the build()
      * function. The placeAtom function does not actually place the atom 
      * in a cell, but calculates and retains a record of the cell index.
      * Computing indices for all atoms before placing them in cells
      * allows the class to compute the amount of memory required for 
      * each cell. The stored cell index for each atom is used to build 
      * the cell list in the build() method. 
      *
      * The method quietly does nothing if the atom is outside the expanded
      * domain for nonbonded ghosts, which extends one cutoff length beyond
      * the domain boundaries the domain boundaries in each direction.
      *
      * Implementation: This function stores the cell index and a pointer 
      * to the Atom in the appropriate element of the tags_ array, and 
      * calls Cell:incrementCapacity() to increment a counter for the
      * number of atoms in the associated Cell. 
      *
      * \param atom  Atom object to be added.
      */
      void placeAtom(Atom &atom);

      /**
      * Build the cell list.
      *
      * This function must be called after completing a loop over atoms,
      * in which placeAtom() is called for each atom. The build() function
      * uses information gathered in this loop to build and fill all of
      * the cells. 
      *
      * Implementation details: This function adds every atom to a Cell,
      * but does not call the CellAtom::update() function to set values 
      * for the position and id. The CellAtom::update() function is 
      * instead later called by CellList::update(). This allows the
      * CellList::build() to be called when atomic coordinates are in
      * scaled [0,1] form, and update() to be called after the positions
      * positions are transformed back to Cartesian coordinates.
      */
      void build();

      /**
      * Update the cell list.
      *
      * This method must be called after the cell list is built, and after
      * transformation back to Cartesian coordinates.
      */
      void update();

      /**
      * Reset the cell list to its empty state (no Atoms).
      */
      void clear();

      /**
      * Return Grid object by const reference.
      */
      const Grid& grid() const;

      /**
      * Dimension of each cell in direction i.
      *
      * \param i Cartesian index i = 0, 1, 2
      */
      double cellLength(int i) const;

      /**
      * Return number of cells that span the cutoff length.
      */
      int nCellCut() const;

      /**
      * Return the index of the cell that contains a position Vector. 
      *
      * Returns a null value of -1 if the position is outside the 
      * expanded domain for nonbonded ghost atoms. 
      *
      * \param position position Vector, inside the boundary
      */
      int cellIndexFromPosition(const Vector& position) const;

      /**
      * Return pointer to first local cell in linked list.
      */
      const Cell* begin() const;

      /**
      * Return a specified cell by const reference.
      * 
      * \param i cell index
      */
      const Cell& cell(int i) const;

      /**
      * Get total number of atoms (local and ghost) in this CellList.
      */
      int nAtom() const;

      /**
      * Get number of atoms that were rejected (not placed in cells)
      */
      int nReject() const;

      #ifdef UTIL_DEBUG
      /**
      * Get maximum number of atoms in one cell.
      */
      int maxNAtomCell() const;
      #endif

      /**
      * Maximum number of atoms for which space is allocated.
      */
      int atomCapacity() const;

      /**
      * Number of cells for which space has been allocated.
      */
      int cellCapacity() const;

      /**
      * Has this CellList been built?
      *
      * IsBuilt is set true by the build() method and false by clear().
      */
      bool isBuilt() const;

      /**
      * Return true if valid, or throw Exception.
      *
      * \return true if valid, otherwise throw an Exception.
      */
      bool isValid() const;

   private:

      /*
      * Temporary storage for atom pointers, before copying to cells. 
      */
      struct Tag {
         Atom* ptr;
         int cellRank;
      };

      /// Array of strips of relative offsets to neighboring cells.
      Cell::OffsetArray offsets_;

      /// Grid for cells.
      Grid grid_;

      /// Array of atom tags (dimension atomCapacity_)
      DArray<Tag> tags_;

      /// Array of CellAtom objects, sorted by cell (dimension atomCapacity_).
      DArray<CellAtom> atoms_;

      /// Array of Cell objects.
      GArray<Cell> cells_;

      /// Lower coordinate bounds (local atoms).
      Vector lower_; 

      /// Upper coordinate bounds (local atoms).
      Vector upper_; 

      /// Length of each cell in grid
      Vector cellLengths_; 

      /// Lower bound for nonbonded ghosts.
      Vector lowerOuter_; 

      /// Upper coordinate bound for nonbonded ghosts.
      Vector upperOuter_; 

      /// Pointer to first local cell (to initialize iterator).
      Cell* begin_;

      /// Total number of atoms in cell list.
      int nAtom_;

      /// Number of atoms that were not placed in cells.
      int nReject_;

      /// Number of cells to span cutoff length
      int nCellCut_;

      #ifdef UTIL_DEBUG
      /// Maximum number of atoms in one cell. 
      int maxNAtomCell_;
      #endif

      /// Has this CellList been built?
      bool isBuilt_;

   }; 

   // Public inline member function definitions

   /*
   * Identify the cell for an Atom, based on its position.
   */
   inline int CellList::cellIndexFromPosition(const Vector& position) const
   {
      IntVector r;
      for (int i = 0; i < Dimension; ++i) {
         if (position[i] <= lowerOuter_[i]) {
            return -1;
         }
         if (position[i] >= upperOuter_[i]) {
            return -1;
         }
         r[i] = int( (position[i] - lowerOuter_[i])/ cellLengths_[i] );
         assert(r[i] < grid_.dimension(i));
      }
      return grid_.rank(r);
   }

   /*
   * Compute atomic cell index, append pointer and index to tags_ array.
   */
   inline void CellList::placeAtom(Atom &atom)
   {
      // Preconditon
      assert(nAtom_ < tags_.capacity());

      int rank = cellIndexFromPosition(atom.position());
      if (rank >= 0) {
         tags_[nAtom_].cellRank = rank;
         tags_[nAtom_].ptr = &atom;
         cells_[rank].incrementCapacity();
         ++nAtom_;
      } else {
         ++nReject_;
      }
   }

   /*
   * Return associated Grid object.
   */
   inline const Grid& CellList::grid() const
   {  return grid_; }

   /*
   * Return length of each cell in direction i.
   */
   inline double CellList::cellLength(int i) const
   {  return cellLengths_[i]; }

   /*
   * Return reference to cell number i.
   */
   inline const Cell& CellList::cell(int i) const
   {
      assert(i < cells_.size());  
      return cells_[i]; 
   }

   /*
   * Return pointer to first Cell.
   */
   inline const Cell* CellList::begin() const
   {  return begin_; }

   /*
   * Get the maximum number of atoms for which space has been allocated.
   */
   inline int CellList::atomCapacity() const
   {  return atoms_.capacity(); }

   /*
   * Get the number of cells for which space has been allocated.
   */
   inline int CellList::cellCapacity() const
   {  return cells_.capacity(); }

   /*
   * Is this CellList built (i.e., full of atoms)?
   */
   inline bool CellList::isBuilt() const
   {  return isBuilt_; }

}
#endif
