#ifndef MCMD_RADIUS_GYRATION_SQ_H
#define MCMD_RADIUS_GYRATION_SQ_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/AverageAnalyzer.h> // base class template
#include <mcMd/simulation/System.h>              // class templ. param.
#include <util/containers/DArray.h>              // member template
#include <util/space/Vector.h>                   // member templ. parameter

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
   * This class calculates the squared radius of gyration by accumulating an
   * average value for the above quantity, and optionally outputs block
   * averages to file at an interval specified by the input parameter
   * nSamplePerBlock. No block averages are output if nSamplePerBlock = 0.
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class RadiusGyrationSq : public AverageAnalyzer<System>
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      RadiusGyrationSq(System &system);
   
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
      * Evaluate squared radii of gyration for all molecules, add to ensemble.
      */
      virtual void compute();
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

   private:

      /// Pointer to array of positions for all beads in a molecule
      DArray<Vector> positions_;
   
      /// Pointer to relevant Species.
      Species* speciesPtr_;
   
      /// Index for relevant species.
      int speciesId_;
   
      /// Number of atoms per molecule of this species.
      int nAtom_;
   
      /// Has readParam been called?
      bool isInitialized_;

   };

   /**
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void RadiusGyrationSq::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & nAtom_;
      ar & positions_;
   }

}
#endif
