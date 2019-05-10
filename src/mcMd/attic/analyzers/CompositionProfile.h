#ifndef MCMD_COMPOSITION_PROFILE_H
#define MCMD_COMPOSITION_PROFILE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>   // base class template
#include <mcMd/simulation/System.h>               // base class template parameter
#include <util/accumulators/Distribution.h>
#include <util/containers/DArray.h>               // member template

#include <util/global.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   /**
   * Distribution of monomer positions along one direction.
   *
   * A CompositionProfile evaluates the distribution of values for
   * one component of the particle position vector. A direction
   * vector is specified as an IntVector containing integer
   * Miller indices.

   * The dot product of monomer position vector and unit
   * direction vector is added to distribution function of
   * particular monomer type and direction vector.
   */
   class CompositionProfile : public SystemAnalyzer<System>
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      CompositionProfile(System &system);

      /**	
      * Destructor.
      */
      ~CompositionProfile();

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int               interval        sampling interval
      *   - string            outputFileName  output file base name
      *   - int               nDirection      number of directions
      *   - DArray<IntVector> intVectors      IntVector directions
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Clear accumulators.
      */
      virtual void setup();

      /**
      * Add particle positions to histogram.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Serialize internal state to/from archive.
      *
      * \param ar       archive
      * \param version  index for archive version
      */
      template <class Archive>
      void serialize(Archive &ar, const unsigned int version);

      /**
      * Output results to predefined output file.
      */
      virtual void output();

   protected:

      /// Output file stream.
      std::ofstream outputFile_;

      /// Distribution statistical accumulators.
      DArray<Distribution> accumulators_;

      /// Array of Miller index vectors for directions.
      DArray<IntVector> intVectors_;

      /// Array of direction vectors.
      DArray<Vector> waveVectors_;

      // Array of log file descriptors
      DArray<std::ofstream> logFiles_;

      /// Number of directions.
      int  nDirection_;

      /// Number of samples thus far.
      int  nSample_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int  nAtomType_;

      /// Number of bins for density profile
      int nBins_;

      /// True if this is the first step
      bool isFirstStep_;

      /// Has readParam been called?
      bool isInitialized_;

      /// Distribution statistical accumulators.
      DArray<Distribution> currentAccumulators_;

      /**
      * Update wavevectors.
      */
      void makeWaveVectors();

   };

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void CompositionProfile::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & nDirection_;
      ar & intVectors_;
      ar & waveVectors_;
      ar & accumulators_;
      ar & nSample_;
      ar & nAtomType_;
   }

}
#endif
