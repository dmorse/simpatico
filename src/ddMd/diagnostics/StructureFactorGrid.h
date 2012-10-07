#ifndef DDMD_STRUCTURE_FACTOR_GRID_H
#define DDMD_STRUCTURE_FACTOR_GRID_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "StructureFactor.h"
#include <util/crystal/LatticeSystem.h>

namespace DdMd
{

   using namespace Util;

   /**
   * StructureFactorGrid evaluates structure factors in Fourier space.
   *
   * This class evaluates the structures factors for all wavevectors
   * for all wavevectors within a grid, within a region in which all
   * of the h, k, l integer wavevector components (i.e., Miller indices) 
   * has an absolute magnitude less than or equal to a parameter hMax.
   * 
   * The class also requires a user to specify a variable lattice of
   * enum type LatticeSystem whose value (e.g., Cubic, Orthorhombic, 
   * Triclinic) determines the symmetry group that is used to group
   * group equivalent wavevectors into "stars". Use Triclinic to 
   * impose no symmetry. 
   *
   * Example: Here is parameter file for a system with nAtomType = 2, 
   * with two monomer types 0 and 1, for calculating the total, 00,
   * and 01 correlation functions on a grid with Miller indices up
   * and including 5:
   *
   * \code
   * StructureFactorGrid{
   *    interval                     1000
   *    outputFileName     StructureFactorGrid
   *    nAtomTypeIdPair                  3
   *    atomTypeIdPairs          -1     -1
   *                              0      0
   *                              0      1
   *    hMax                             5
        lattice                      Cubic
   * }
   * \endcode
   * At the end of a simulation, all of the structure factors are
   * output in a file with a suffix *.dat. Each line in this file
   * contains the 3 Miller indices of a wavevector, the absolute
   * magnitude of the wavevector, and a list of nAtomTypeIdPair 
   * structure factor values for the wavevector, one for each 
   * atomTypeId pair.
   * 
   * \ingroup DdMd_Diagnostic_Module
   */
   class StructureFactorGrid : public StructureFactor
   {

   public:

      /**	
      * Constructor.
      *
      * \param simulation reference to parent DdMd::Simulation object
      */
      StructureFactorGrid(Simulation &simulation);

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int               interval        sampling interval 
      *   - string            outputFileName  output file base name
      *   - int               nMode           number of modes
      *   - DMatrix<double>   modes           mode vectors
      *   - int               hMax            maximum Miller index
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Output structure factors, averaged over stars.
      */
      virtual void output();

   private:

      /// Array of ids for first wavevector in each star.
      DArray<int>  starIds_;

      /// Array of star sizes.
      DArray<int>  starSizes_;

      /// Maximum Miller index of wavevectors in grid.
      int   hMax_;

      /// Number of stars of symmetry related wavevectors.
      int   nStar_;

      /// Lattice system used to create stars.
      LatticeSystem   lattice_;

   };

}
#endif
