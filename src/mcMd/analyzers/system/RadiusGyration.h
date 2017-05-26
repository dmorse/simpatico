#ifndef MCMD_RADIUS_GYRATION_H
#define MCMD_RADIUS_GYRATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
#include <mcMd/simulation/System.h>             // class template parameter
#include <util/accumulators/Average.h>          // member
#include <util/containers/DArray.h>             // member template
#include <util/space/Vector.h>                  // member template parameter

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
   * Radius of gyration of a molecule.
   *
   * The radius of gyration R_g is the root mean-squared separation of 
   * atoms in a molecule from the molecular center of mass.  For a species
   * of molecule with N atoms, with positions \f$ R_0, R_1, ...., R_{N-1} \f$,
   * we define
   * \f[
   *  
   *   R_g^2 = \sum_i \langle | R_i - R_{cm} |^2 \rangle / N
   * 
   * \f]
   * where \f$ R_{cm} = (R_0 + R_1 + ... R_{N-1})/N \f$ is the center of mass, 
   * and where \f$ \langle \cdots \rangle \f$ denotes an ensemble average over 
   * all molecules of a specified species. 
   *
   * This class calculates the radius of gyration by accumulating an
   * average value for the above quantity, and optionally outputs block
   * averages to file at an interval specified by the input parameter
   * nSamplePerBlock. No block averages are output if nSamplePerBlock = 0.
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class RadiusGyration : public SystemAnalyzer<System>
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      RadiusGyration(System &system);
   
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
      * Load state to an archive.
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

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Average object - statistical accumulator
      Average  accumulator_;

      /// Pointer to array of positions for all beads in a molecule
      DArray<Vector> positions_;
   
      /// Pointer to relevant Species.
      Species *speciesPtr_;
   
      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Index for relevant species.
      int     speciesId_;
   
      /// Number of atoms per molecule of this species.
      int     nAtom_;
   
      /// Has readParam been called?
      bool    isInitialized_;

   };

   /**
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void RadiusGyration::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & nSamplePerBlock_;
      ar & speciesId_;
      ar & nAtom_;
      ar & accumulator_;
      ar & positions_;
   }

}
#endif
