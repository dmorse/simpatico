#ifndef MCMD_INTRA_STRUCTURE_FACTOR_GRID_H
#define MCMD_INTRA_STRUCTURE_FACTOR_GRID_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraStructureFactor.h"
#include <util/crystal/LatticeSystem.h>

namespace McMd
{

   using namespace Util;

   /**
   * IntraStructureFactorGrid evaluates structure factors in Fourier space.
   *
   * This class evaluates the intramolecular contribution to the 
   * structures factor for all wavevectors within a grid,
   * within a region in which all of the h, k, l integer wavevector
   * components (i.e., Miller indices) has an absolute magnitude less than
   * or equal to a parameter hMax.
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
   * and including 5 for molecules of species 0:
   *
   * \code
   * IntraStructureFactorGrid{
   *    interval                     1000
   *    outputFileName     IntraStructureFactorGrid
   *    speciesId                        0
   *    nAtomTypeIdPair                            3
   *    atomTypeIdPairs                    -1     -1
   *                              0      0
   *                              0      1
   *    hMax                             5
   *    lattice                      Cubic
   * }
   * \endcode
   * At the end of a simulation, all of the structure factors are
   * output in a file with a suffix _avg.dat. Each line in this file
   * contains the 3 Miller indices of a wavevector, the absolute
   * magnitude of the wavevector, and a list of nAtomTypeIdPair 
   * structure factor values for the wavevector, one for each 
   * atomTypeId pair. During the simulation, at each sampling, a log file
   * with suffix .dat containing the instantaneous values of the
   * structure factors is appended to.
   * 
   * \ingroup McMd_Diagnostic_Module
   */
   class IntraStructureFactorGrid : public IntraStructureFactor
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      IntraStructureFactorGrid(System &system);

      /**	
      * Destructor.
      */
      ~IntraStructureFactorGrid();

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int                 interval         sampling interval 
      *   - int                 speciesId        species identifier
      *   - string              outputFileName   output file base name
      *   - int                 nAtomTypeIdPair  number of atomTypeIdPairs
      *   - DArray< Pair<int> > atomTypeIdPairs  atomTypeIdPair vectors
      *   - int                 hMax             maximum Miller index
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Set up before a simulation.
      */
      virtual void setup();

      /**
      * Output structure factors, averaged over stars.
      */
      virtual void output();

      /**
      * Add particles to StructureFactor accumulators.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);

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

      /// Log file
      std::ofstream logFile_;

      /// Is this the first step?
      bool isFirstStep_;

      /// Has readParam been called?
      bool isInitialized_;
   };


   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void IntraStructureFactorGrid::serialize(Archive& ar, const unsigned int version)
   {
      IntraStructureFactor::serialize(ar, version);
      ar & hMax_;
      //serializeEnum(ar, lattice_);
      ar & lattice_;
      ar & nStar_;
      ar & starIds_;
      ar & starSizes_;
   }

}
#endif
