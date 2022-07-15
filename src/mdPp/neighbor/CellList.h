#ifndef MDPP_CELL_LIST_H
#define MDPP_CELL_LIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mdPp/neighbor/Cell.h>
#include <mdPp/neighbor/CellAtom.h>
#include <mdPp/chemistry/Atom.h>
#include <simp/boundary/Boundary.h>
#include <util/space/Grid.h>
#include <util/containers/DArray.h>
#include <util/containers/GArray.h>
#include <util/containers/FArray.h>
#include <util/global.h>

namespace MdPp
{

   using namespace Util;
   using namespace Simp;

   /**
   * A cell list used only to identify nearby atom pairs.
   *
   * An CellList divides the simulation unit cell into a grid of cells,
   * such that the length of each cell in each direction is greater than 
   * a specified cutoff distance. 
   *
   * Building an CellList (Cartesian coordinates);
   * \code
   *
   *    AtomStorage storage;
   *    CellList cellList;
   *    Vector lower;        // Vector of lower bounds (local atoms)
   *    Vector upper;        // Vector of upper bounds (local atoms)
   *    Vector cutoffs;      // Vector of cutoff lengths for each axis.
   *    int atomCapacity     // max number of atoms on this processor
   *
   *    // Set elements of cutoffs vector to same value
   *    for (int i = 0; i < Dimension; ++i) {
   *       cutoffs[i] = cutoff;
   *    } 
   *
   *    // Bounds on lower and upper used here to allocate memory.
   *    cellList.allocate(atomCapacity, lower, upper, cutoffs);
   *  
   *    // Make the actual grid and clear it.
   *    cellList.makeGrid(lower, upper, cutoffs);
   *    cellList.clear();
   *
   *    // Place all atoms.
   *    AtomStorage::AtomIterator  atomIter;
   *    for (storage.begin(atomIter); atomIter.notEnd(); ++atomIter) {
   *       cellList.placeAtom(*atomIter);
   *    }
   *
   *    // Build cell list
   *    cellList.build();
   *
   * \endcode
   *
   * The atomCapacity parameter should be set equal to the sum of atomCapacity 
   * and the ghostCapacity of the associated atomStorage, which is the maximum 
   * total number of atoms that can exist on this processor.
   *
   * If the upper and lower bounds and atom coordinates are all expressed in
   * generalized coordinates, which span 0.0 - 1.0 over the primitive periodic
   * cell in each direction, each element of the cutoffs vector is given by a
   * ratio cutoffs[i] = cutoff/length[i], where length[i] is the Cartesian
   * distance across the unit cell along a direction parallel to reciprocal 
   * basis vector i. 
   *
   * See Cell documentation for an example of how to iterate over local cells 
   * and neighboring atom pairs. 
   * 
   * \ingroup MdPp_Neighbor_Module
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
      */
      virtual ~CellList();

      /**
      * Allocate memory for this CellList (generalized coordinates).
      *
      * This function:
      *
      *   - Allocates an array of atomCapacity CellList::Tag objects.
      *   - Allocates an array of atomCapacity CellAtom objects.
      *   - Allocates an array of Cell objects sized for this boundary.
      *
      * The elements of the lower, upper, and cutoffs parameters should 
      * contain the lower and upper coordinate bounds for this processor, 
      * and cutoff values in each direction, for a boundary was chosen to
      * be larger than any that will be encountered during the simulation
      * These parameters are used only to allocate memory.
      *
      * This version of the function is designed for use with generalized
      * coordinates. See the makeGrid() method for a discussion of the
      * parameter values in generalized coordinates.
      *
      * \param atomCapacity dimension of global array of atoms
      * \param lower        lower coordinate bounds for this processor
      * \param upper        upper coordinate bounds for this processor
      * \param cutoffs      pair cutoff distance in each direction
      */
      void allocate(int atomCapacity, const Vector& lower, const Vector& upper, 
                    const Vector& cutoffs);

      /**
      * Allocate memory for this CellList (Cartesian coordinates).
      *
      * This function is designed for use with Cartesian coordinates, for which
      * the lower, upper and cutoff parameters all have dimensions of length.
      * The function calls the allocate() method with a Vector of cutoffs
      * internally, after setting every element of the Vector to the same value.
      *
      * \param atomCapacity dimension of global array of atoms
      * \param lower        lower bound for this processor in maximum boundary
      * \param upper        upper bound for this processor in maximum boundary
      * \param cutoff       pair cutoff distance in each direction
      */
      void allocate(int atomCapacity, const Vector& lower, const Vector& upper, 
                    double cutoff);

      /**
      * Make the cell grid (using generalized coordinates).
      *
      * This method makes a cell grid in which the number of cells in each
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
      * distance across the primitive unit cell along the direction parallel 
      * to reciprocal lattice basis vector i.
      *
      * \param lower    lower bound of local atom coordinates.
      * \param upper    upper bound of local atom coordinates.
      * \param cutoffs  pair cutoff length in each direction
      */
      void 
      makeGrid(const Vector& lower, const Vector& upper, const Vector& cutoffs);

      /**
      * Determine the appropriate cell for an Atom, based on its position.
      *
      * This method does not place the atom in a cell, but calculates a 
      * cell index and retains the value, which is used to place atoms in
      * the build() method.
      *
      * The method quietly does nothing if the atom is outside the domain.
      *
      * \param atom  Atom object to be added.
      */
      void placeAtom(Atom &atom);

      /**
      * Build the cell list.
      *
      * This method must be called after completing a loop over atoms,
      * in which placeAtom() is called for each atom. The build() method
      * uses information gathered in this loop to build and fill all of
      * the cells. 
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
      * Reset the cell list to its empty state (no atoms).
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
      * Has memory been allocated for this CellList?
      */
      bool isAllocated() const;

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

      /// Pointer to first local cell (to initialize iterator).
      Cell* begin_;

      /// Total number of atoms in cell list.
      int nAtom_;

      /// Number of atoms that were not placed in cells.
      int nReject_;

      #ifdef UTIL_DEBUG
      /// Maximum number of atoms in one cell. 
      int maxNAtomCell_;
      #endif

      /// Has this CellList been built?
      bool isBuilt_;

      /// Stores the same offsets as the grid. Used to calculate Cell::OffsetArray.
      int gridOffsets_ [3];

      /**
      * Calculate required dimensions for cell grid and resize cells_ array.
      *
      * Called internally by allocate and makeGrid. Resizes cells_ array to 
      * match size of new grid. Does not link cells or calculate offsets to 
      * neighbors.
      *
      * \param lower  lower bound used to allocate array of cells.
      * \param uppper  upper bound used to allocate array of cells.
      * \param cutoffs  minimum dimension of a cell in each direction
      */
      void setGridDimensions(const Vector& lower, const Vector& upper,
                             const Vector& cutoffs);

      /**
      * Return true if atomId is valid, i.e., if 0 <= 0 < atomCapacity.
      */
      bool isValidAtomId(int atomId);

      /**
      * Return the offset such that (&cells_[cellId] + offset) is a pointer
      * to the cell a distance x from cells_[cellId] along the axis i shifted 
      * within the range 0 <= x < grid_.Dimension(i)
      *
      * \param cellId id of cell to calculate offset from
      * \param i index for Cartesian axis. Must be 0, 1, or 2.
      * \param x integer coordinate of cell along axis i.
      */
      int calculateAxisOffset(int cellId, int i, int x) const;
   };

   // Private inline method definitions:

   /*
   * Return value of offset for cell a distance x along axis i shifted to the range
   * 0 <= x <= grid_.Dimension(i)
   */
   inline int CellList::calculateAxisOffset(int cellId, int i, int x) const
   {
      // mod the shift so it's within the grid dimension
      x %= grid_.dimension(i);

      // calculate cell's position within the grid
      int p[3]; int j;
      for (j = 0; j < Dimension -1; ++j) {
         p[j] = cellId/gridOffsets_[j];
         cellId -= p[j]*gridOffsets_[j];
      }
      p[j] = cellId;

      if (p[i] + x >= grid_.dimension(i)) return (x - grid_.dimension(i))*gridOffsets_[i];
      if (p[i] + x < 0) return (grid_.dimension(i) + x)*gridOffsets_[i];
      return x*gridOffsets_[i];
   }

   // Public inline method definitions:

   /*
   * Identify the cell within which a position lies.
   */
   inline int CellList::cellIndexFromPosition(const Vector& position) const
   {
      IntVector r;
      for (int i = 0; i < Dimension; ++i) {
         if (position[i] <= lower_[i]) {
            return -1;
         }
         if (position[i] >= upper_[i]) {
            return -1;
         }
         r[i] = int( (position[i] - lower_[i]) / cellLengths_[i] );
         assert(r[i] < grid_.dimension(i));
      }
      return grid_.rank(r);
   }

   /*
   * Add an atom to the appropriate cell, based on its position.
   */
   inline void CellList::placeAtom(Atom &atom)
   {
      // Preconditon
      assert(nAtom_ < tags_.capacity());

      int rank = cellIndexFromPosition(atom.position);
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
   * Return true iff atomId is valid, i.e., if 0 <= 0 < atomCapacity.
   */
   inline bool CellList::isValidAtomId(int atomId)
   {  return ( (0 <= atomId) && (atomId < tags_.capacity()) ); }

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
   {  return cells_[i]; }

   /*
   * Return pointer to first Cell.
   */
   inline const Cell* CellList::begin() const
   {  return begin_; }

   /*
   * Is this CellList allocated?
   */
   inline bool CellList::isAllocated() const
   {  return (cells_.capacity() > 0); }

}
#endif
