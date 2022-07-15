#ifndef MCMD_CLUSTER_HISTOGRAM_H
#define MCMD_CLUSTER_HISTOGRAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>           // base class templ
#include <mcMd/simulation/System.h>                  // class templ param
#include <mcMd/analyzers/system/ClusterIdentifier.h> // member
#include <util/accumulators/IntDistribution.h>       // member

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * Identify micelle clusters in polymeric systems.
   */
   class ClusterHistogram : public SystemAnalyzer<System>
   {

   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      ClusterHistogram(System &system);

      /**
      * Read parameters from file, and allocate data array.
      *
      * Input format:
      *
      *   - int    interval         : sampling interval
      *   - string outputFileName   : base name for output file(s)
      *   - int    speciesId        : integer id for Species of interest
      *   - int    atomTypeId       : integer id for core atom type
      *   - double cutoff           : distance cutoff
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Clear accumulator.
      */
      virtual void setup();

      /**
      * Identify clusters in configuration.
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

      /// ClusterIdentifier
      ClusterIdentifier identifier_;

      /// Distribution of the Clusters.
      IntDistribution  hist_;

      /// Output file stream
      std::ofstream outputFile_;

      /// Clusters species type
      int  speciesId_;

      /// Clusters species type
      int  atomTypeId_;

      /// Distance cutoff
      double cutoff_;

      /// Histogram minimum value
      int  histMin_;

      /// Histogram maximum value.
      int  histMax_;

      /// Number of configurations dumped thus far (first dump is zero).
      long  nSample_;

      /// Has readParam been called?
      bool  isInitialized_;

   };

   /**
   * Serialize to/from an archive.
   */
   template <class Archive>
   void ClusterHistogram::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & histMin_;
      ar & histMax_;
      ar & nSample_;
      ar & hist_;
   }

}
#endif
