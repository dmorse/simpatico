#ifndef MCMD_BLOCK_RADIUS_GYRATION_H
#define MCMD_BLOCK_RADIUS_GYRATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h> // base class template
#include <mcMd/simulation/System.h>        // class template parameter
#include <util/accumulators/Average.h>     // member
#include <util/containers/DArray.h>        // member template
#include <util/space/Vector.h>             // member template parameter

#include <cstdio> 
#include <cstring> 

namespace Simp {
   class Species;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * Radius of gyration of different blocks in a molecule.
   *
   * The radius of gyration R_{g,i} of a block of atom type i is the 
   * root mean-squared separation of atoms in a block from the center 
   * of mass of block.  For a species of molecule with N_i atoms of 
   * type i, with positions \f$ R_0, R_1, ...., R_{N_{i}-1} \f$, we 
   * define
   * \f[
   *  
   *   R_{g,i}^2 = \sum_j \langle | R_j - R_{cm,i} |^2 \rangle / N_i
   * 
   * \f]
   * where \f$ R_{cm,i} = (R_0 + R_1 + ... R_{N_{i}-1})/N_i \f$ is the center of mass 
   * of block, and where \f$ \langle \cdots \rangle \f$ denotes an ensemble average 
   * over all molecules of a specified species. 
   *
   * This class calculates the radius of gyration of blocks by accumulating an
   * average value for the above quantity, and optionally outputs block
   * averages to file at an interval specified by the input parameter
   * nSamplePerBlock. No block averages are output if nSamplePerBlock = 0.
   *
   * \sa \ref mcMd_analyzer_BlockRadiusGyration_page "parameter file format"
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class BlockRadiusGyration : public SystemAnalyzer<System>
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      BlockRadiusGyration(System &system);
   
      /**
      * Read parameters from file, and allocate data array.
      *
      * Input format:
      *
      *   - int    interval        : sampling interval
      *   - string outputFileName  : base name for output file(s)
      *   - int    nSamplePerBlock : interval for output of block averages
      *   - int    speciesId       : integer id for Species of interest
      *
      * No block averages are output if nSamplePerBlock = 0. Otherwise,
      * block averages are output to a file named (outputFileName).dat. 
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);
   
      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /**
      * Evaluate squared radii of gyration for all molecules, add to ensemble.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Serialize to/from an archive. 
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Array of Average objects - statistical accumulators.
      DArray<Average>  accumulators_;

      /// Array of positions for all atoms in a molecule
      DArray<Vector> positions_;
      
      /// Array of center of mass vectors for blocks of different types
      DArray<Vector> rCom_;
      
      /// Array of number of atoms in blocks of different types
      DArray<int>  iTypeNAtom_;
      
      /// Values of dRSq for atom types i=j (temporary storage)
      DArray<double> dRSq_; 

      /// Values of dRSq for atom types i !=j (not persistent)
      DArray<double> dRSqPair_;  

      /// Pointer to relevant Species.
      Species*  speciesPtr_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int  nAtomType_;
      
      /// Number of distinct types of type pairs (this is 1 in a diblock melt)
      int  nAtomTypePairs_;
      
      /// Number of samples per block average output.
      int  nSamplePerBlock_;

      /// Pointer to relevant Species.
      int  speciesId_;
   
      /// Number of atoms per molecule of this species.
      int  nAtom_;
   
      /// Has readParam been called?
      bool  isInitialized_;

   };

   /**
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void BlockRadiusGyration::serialize(Archive& ar, const unsigned int version)
   {  
      Analyzer::serialize(ar, version);
      ar & nSamplePerBlock_;
      ar & speciesId_;
      ar & nAtom_;
      ar & nAtomType_;
      ar & nAtomTypePairs_;
      ar & accumulators_;
   }

}
#endif
